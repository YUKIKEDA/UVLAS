param([string]$Scratch, [string]$Repo, [string]$Packages)
$ErrorActionPreference = 'Stop'
$build = Join-Path $Repo 'out/issue166-build'
$source = Join-Path $Repo 'src/aero/pchip_slot_integer.cpp'
$raw = [IO.File]::ReadAllBytes($source)
$original = [Text.Encoding]::UTF8.GetString($raw).Replace("`r`n","`n")
$sourceHash = [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($raw))
if ($sourceHash -ne '7CE3807A949B725E003D074EE5FE170A8642856353798EAE357E43E2EAA7BBD5') { throw 'Unexpected mutation source' }
function Set-SingleLiteral([string]$text, [string]$before, [string]$after) {
    if ([regex]::Matches($text,[regex]::Escape($before)).Count -ne 1) { throw "Mutation site is not unique: $before" }
    $text.Replace($before,$after)
}
$noFloor = Set-SingleLiteral $original 'if (in.sign < 0 && r.used() != 0) {' 'if (in.sign > 0 && r.used() != 0) {'
$carryLine = '            out.words[i + j] = static_cast<std::uint32_t>(v);' + [char]10 + '            carry = v >> 32U;'
$noCarry = Set-SingleLiteral $original $carryLine ('            out.words[i + j] = static_cast<std::uint32_t>(v);' + [char]10 + '            carry = 0;')
$inputCases = @(
    @{name='control'; source=$original; filter='Private PCHIP slots preserve exact signed rational arithmetic,Private PCHIP dense products and cross products match independent big integers'},
    @{name='wrong_floor'; source=$noFloor; filter='Private PCHIP slots preserve exact signed rational arithmetic'},
    @{name='lost_carry'; source=$noCarry; filter='Private PCHIP dense products and cross products match independent big integers'}
)
Write-Output "source_sha256=$sourceHash"
foreach($case in $inputCases) {
    $src = Join-Path $Scratch ($case.name + '.cpp')
    $obj = Join-Path $Scratch ($case.name + '.obj')
    $exe = Join-Path $Scratch ($case.name + '.exe')
    [IO.File]::WriteAllBytes($src,[Text.Encoding]::UTF8.GetBytes($case.source))
    $compile = @('/nologo','/c','/MDd','/Od','/EHsc','/std:c++23preview','/permissive-','/utf-8','/W4','/WX','/DWIN32_LEAN_AND_MEAN','/DNOMINMAX',"/I$Repo/include","/I$Repo/src/aero","/external:I$Packages/include",'/external:W0',"/Fo$obj",$src)
    & cl @compile
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $test = Join-Path $build 'tests/CMakeFiles/uvlas_test_aero_pchip_rational_slots.dir/Debug/aero/pchip_rational_slots_test.cpp.obj'
    $linkArgs = @('/nologo','/subsystem:console','/machine:x64',"/out:$exe",$obj,$test,"$Packages/debug/lib/manual-link/Catch2Maind.lib","$Packages/debug/lib/Catch2d.lib","$build/Debug/uvlas_aero.lib","$build/Debug/uvlas_core.lib",'psapi.lib')
    & link @linkArgs
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $check = & $exe $case.filter --reporter compact
    $code = $LASTEXITCODE
    if ($case.name -eq 'control') {
        if ($code -ne 0 -or ($check -join "`n") -notmatch 'passed') { throw 'Unmutated control did not pass' }
    } elseif ($code -eq 0 -or ($check -join "`n") -notmatch 'failed:') {
        $check | Select-Object -First 5
        throw "Mutation did not trigger the expected assertion: $($case.name)"
    }
    $relevant = @($check | Where-Object { $_ -match 'assertions?|passed' })
    Write-Output ("{0}: exit={1}; {2}" -f $case.name,$code,($relevant -join ' '))
    if ($code -ne 0) { $check | Where-Object {$_ -match 'failed:' } | Select-Object -First 1 | ForEach-Object { if ($_ -match 'failed: (.{0,150})') { $Matches[1] } } }
}
