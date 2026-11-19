param([Parameter(Mandatory)][string]$Repo,[Parameter(Mandatory)][string]$Build,[Parameter(Mandatory)][string]$Packages,[Parameter(Mandatory)][string]$Scratch)
$ErrorActionPreference='Stop'
# Reproduction only: use the frozen numeric source and matching Debug objects at b16c908.
$revision=(& git -C $Repo rev-parse HEAD).Trim()
if($LASTEXITCODE -ne 0 -or $revision -ne 'b16c908f9ea177952f615f462bbc876f20dfebcf'){throw 'Use the frozen b16c908 source tree and matching build'}
if(Test-Path -LiteralPath $Scratch){if(@(Get-ChildItem -LiteralPath $Scratch).Count -ne 0){throw 'Scratch must be empty'}}else{New-Item -ItemType Directory -Path $Scratch|Out-Null}
$raw=[IO.File]::ReadAllBytes((Join-Path $Repo 'src/aero/pchip_coefficients.cpp'))
$sha=[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($raw))
if($sha -ne '843E4D5FA5B2D6F6DD8BF48AA66FA3D7A6718D92F4CD2BD80F55E46DE07F83FB'){throw 'Unexpected numeric source bytes'}
$original=[Text.Encoding]::UTF8.GetString($raw)
function Replace-SingleSite([string]$source,[string]$before,[string]$after){if([regex]::Matches($source,[regex]::Escape($before)).Count -ne 1){throw "Mutation site is not unique: $before"};return $source.Replace($before,$after)}
$cases=@(
    @{name='control';source=$original},
    @{name='cubic_sign';source=(Replace-SingleSite $original 'm.multiply(m.integer(2), m.subtract(a, b))' 'm.multiply(m.integer(2), m.subtract(b, a))')},
    @{name='endpoint_zero';source=(Replace-SingleSite $original 'sz != sd,' 'sz == sd,')}
)
foreach($case in $cases){$source=Join-Path $Scratch ($case.name+'.cpp');$object=Join-Path $Scratch ($case.name+'.obj');$exe=Join-Path $Scratch ($case.name+'.exe')
    [IO.File]::WriteAllBytes($source,[Text.Encoding]::UTF8.GetBytes($case.source))
    & cl /nologo /c /MDd /Od /EHsc /std:c++23preview /permissive- /utf-8 /W4 /WX /DWIN32_LEAN_AND_MEAN /DNOMINMAX "/I$Repo/include" "/I$Repo/src/aero" "/external:I$Packages/include" /external:W0 "/Fo$object" $source
    if($LASTEXITCODE -ne 0){throw "Compile failed: $($case.name)"}
    $test=Join-Path $Build 'tests/CMakeFiles/uvlas_test_aero_pchip_coefficients.dir/Debug/aero/pchip_coefficients_test.cpp.obj'
    & link /nologo /subsystem:console /machine:x64 "/out:$exe" $object $test "$Packages/debug/lib/manual-link/Catch2Maind.lib" "$Packages/debug/lib/Catch2d.lib" "$Build/Debug/uvlas_aero.lib" "$Build/Debug/uvlas_core.lib" psapi.lib
    if($LASTEXITCODE -ne 0){throw "Link failed: $($case.name)"}
    $output=& $exe 'Private PCHIP coefficients match frozen exact witnesses and independent Hermite equations' --rng-seed 176 --reporter compact
    $code=$LASTEXITCODE
    $expected=if($case.name -eq 'control'){0}else{42}
    if($code -ne $expected -or ($case.name -ne 'control' -and ($output -join "`n") -notmatch 'failed:')){throw "Unexpected test result: $($case.name), exit=$code"}
    $output|Set-Content -LiteralPath (Join-Path $Scratch ($case.name+'.test.txt')) -Encoding utf8
    Write-Output "$($case.name): compile=0 link=0 Catch=$code"
}
exit 0
