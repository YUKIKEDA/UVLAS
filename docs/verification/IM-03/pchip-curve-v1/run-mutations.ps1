param([string] $BuildDirectory = 'out/issue166-build')
$ErrorActionPreference = 'Stop'
$project = $PWD.Path
$build = Join-Path $project $BuildDirectory
$original = Join-Path $project 'src/aero/pchip_curve.cpp'
$vcpkg = ((Select-String -LiteralPath "$build/CMakeCache.txt" -Pattern '^VCPKG_INSTALLED_DIR:PATH=').Line.Split('=',2)[1]) + '/x64-windows'
$scratch = Join-Path ([IO.Path]::GetTempPath()) ("uvlas-258-curve-" + [guid]::NewGuid().ToString('N'))
[void][IO.Directory]::CreateDirectory($scratch)
$generated = [System.Collections.Generic.List[string]]::new()
function New-Artifact($name) {
    $path = Join-Path $scratch $name
    $generated.Add($path)
    return $path
}
try {
    $msvc = Join-Path $project 'tools/dev/msvc.ps1'
    $flags = @('/nologo','/TP','/c','/std:c++23preview','/MDd','/EHsc','/Od','/W4','/WX','/permissive-','/utf-8','/Zc:__cplusplus','/fp:strict','/DNOMINMAX','/DWIN32_LEAN_AND_MEAN',
               "/I$project/src","/I$project/src/aero","/I$project/include","/external:I$vcpkg/include",'/external:W0')
    $test = "$build/tests/CMakeFiles/uvlas_test_aero_pchip_curve.dir/Debug/aero/pchip_curve_test.cpp.obj"
    $libs = @("$vcpkg/debug/lib/manual-link/Catch2Maind.lib","$build/Debug/uvlas_aero.lib","$vcpkg/debug/lib/Catch2d.lib","$build/Debug/uvlas_core.lib",'psapi.lib')
    foreach ($mode in @('control','skip-saved-coefficients','accept-unequal-automatic-sides')) {
        $source = New-Artifact "$mode.cpp"
        [IO.File]::Copy($original,$source)
        if ((Get-FileHash -LiteralPath $source).Hash -ne (Get-FileHash -LiteralPath $original).Hash) { throw 'Product copy differs' }
        if ($mode -ne 'control') {
            $sourceText = [IO.File]::ReadAllText($source)
            if ($mode -eq 'skip-saved-coefficients') {
                $old = 'if (comparison != 0) {'
                $new = 'if (comparison == 0) {'
            } else {
                $old = 'const auto agreement = agree_side(*both, *left, *right, plan.agreement, stop, details);'
                $new = 'static_cast<void>(*both); const core::Result<SideAgreement> agreement = SideAgreement{JetCandidate::Right, 0};'
            }
            if ([regex]::Matches($sourceText,[regex]::Escape($old)).Count -ne 1) { throw ("Unexpected mutation anchor: " + $old) }
            [IO.File]::WriteAllText($source, $sourceText.Replace($old,$new))
        }
        $object = New-Artifact "$mode.obj"
        $exe = New-Artifact "$mode.exe"
        pwsh -NoProfile -File $msvc cl @flags "/Fo$object" $source
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
        pwsh -NoProfile -File $msvc link /nologo /subsystem:console "/out:$exe" $test $object @libs
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
        $savedPath = $env:PATH
        try {
            $env:PATH = "$vcpkg/debug/bin;$savedPath"
            $filter = switch ($mode) {
                'skip-saved-coefficients' { 'PCHIP snapshot retains*' }
                'accept-unequal-automatic-sides' { 'Qualified PCHIP policies*' }
                default { 'PCHIP snapshot retains*,Qualified PCHIP policies*' }
            }
            & $exe $filter --rng-seed 258 --reporter compact
            $outcome = $LASTEXITCODE
        } finally {
            $env:PATH = $savedPath
        }
        $want = if ($mode -eq 'control') { 0 } else { 42 }
        if ($outcome -ne $want) { throw "Unexpected $mode exit $outcome" }
        "$mode Catch exit $outcome"
    }
} finally {
    foreach ($path in $generated) {
        if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path }
    }
    Remove-Item -LiteralPath $scratch
}
