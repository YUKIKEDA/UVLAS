param([string] $BuildDirectory = 'out/issue166-build')
$ErrorActionPreference = 'Stop'
$project = $PWD.Path
$build = Join-Path $project $BuildDirectory
$original = Join-Path $project 'src/aero/p1_shared_action.cpp'
$vcpkg = ((Select-String -LiteralPath "$build/CMakeCache.txt" -Pattern '^VCPKG_INSTALLED_DIR:PATH=').Line.Split('=',2)[1]) + '/x64-windows'
$scratch = Join-Path ([IO.Path]::GetTempPath()) ("uvlas-266-p1-" + [guid]::NewGuid().ToString('N'))
[void][IO.Directory]::CreateDirectory($scratch)
$generated = [System.Collections.Generic.List[string]]::new()
function New-Artifact($name) {
    $path = Join-Path $scratch $name
    $generated.Add($path)
    return $path
}
function Set-Exact($inputText, $pattern, $replacement) {
    if ([regex]::Matches($inputText,$pattern).Count -ne 1) { throw "Unexpected mutation anchor $pattern" }
    return [regex]::Replace($inputText,$pattern,$replacement)
}
try {
    $msvc = Join-Path $project 'tools/dev/msvc.ps1'
    $flags = @('/nologo','/TP','/c','/std:c++23preview','/MDd','/EHsc','/Od','/W4','/WX','/permissive-','/utf-8','/Zc:__cplusplus','/fp:strict','/DNOMINMAX','/DWIN32_LEAN_AND_MEAN',
               "/I$project/src","/I$project/src/aero","/I$project/include","/external:I$vcpkg/include",'/external:W0')
    $test = "$build/tests/CMakeFiles/uvlas_test_aero_p1_shared_action.dir/Debug/aero/p1_shared_action_test.cpp.obj"
    $libs = @("$vcpkg/debug/lib/manual-link/Catch2Maind.lib","$build/Debug/uvlas_aero.lib","$vcpkg/debug/lib/Catch2d.lib","$build/Debug/uvlas_core.lib",'psapi.lib')
    foreach ($mode in @('control','ignore-trace','omit-point-derivative')) {
        $source = New-Artifact "$mode.cpp"
        [IO.File]::Copy($original,$source)
        if ((Get-FileHash -LiteralPath $source).Hash -ne (Get-FileHash -LiteralPath $original).Hash) { throw 'Product copy differs' }
        $sourceText = [IO.File]::ReadAllText($source)
        if ($mode -eq 'ignore-trace') {
            $sourceText = Set-Exact $sourceText 'edge\.owner, edge\.trace, edge\.positive' 'edge.owner, edge.positive'
            $sourceText = Set-Exact $sourceText '\s*&&\s*edges\[last\]\.trace == edges\[first\]\.trace' ''
            [IO.File]::WriteAllText($source, $sourceText)
        }
        if ($mode -eq 'omit-point-derivative') {
            $sourceText = Set-Exact $sourceText '(append\(a, query\[a\]\);\s*for \(std::size_t point = 0; point < )3' '${1}2'
            [IO.File]::WriteAllText($source, $sourceText)
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
            & $exe 'Shared P1 quad action and all derivatives*' --rng-seed 266 --reporter compact
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
