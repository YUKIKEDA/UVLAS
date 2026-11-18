param([string]$Root=(Get-Location).Path,
      [string]$CompileDatabase='out/issue127-tidy-db/compile_commands.json',
      [string]$Build='out/issue127-build',
      [string]$Output='out/p1-release-profile-reproduction')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $Root
$entry=(Get-Content $CompileDatabase -Raw | ConvertFrom-Json)[0]
$include=@($entry.arguments | Where-Object {$_ -like '-external:I*'})[0]
$prefix=Split-Path -Parent $include.Substring('-external:I'.Length)
$source=Join-Path $Root 'docs/verification/IM-03/p1-kernel-implementation-v1'
$destination=[IO.Path]::GetFullPath((Join-Path $Root $Output))
if(-not $destination.StartsWith((Join-Path $Root 'out/'),[StringComparison]::OrdinalIgnoreCase)){
    throw 'Output must stay under the workspace out directory.'
}
if(Test-Path $destination){throw 'Output must be new.'}
foreach($variant in @('before','after')){
    $dir=Join-Path $destination $variant
    New-Item -ItemType Directory -Path "$dir/objects" -Force | Out-Null
    Copy-Item -LiteralPath "$source/profile-$variant.cpp" -Destination "$dir/p1_induction.cpp"
    Copy-Item -LiteralPath "$source/profile-release-$variant.hpp" -Destination "$dir/p1_induction.hpp"
    Copy-Item -LiteralPath "$source/profile-release-main.cpp" -Destination "$dir/profile-main.cpp"
    $arguments=@('/nologo','/std:c++23preview','/EHsc','/O2','/MD','/fp:strict','/utf-8',
        '/DNOMINMAX','/DWIN32_LEAN_AND_MEAN',"/I$dir",'/Isrc','/Isrc/aero','/Iinclude','/Itests',
        "/external:I$prefix/include",'/external:W0',"/Fo$dir/objects/","/Fd$dir/compile.pdb",
        "$dir/p1_induction.cpp","$dir/profile-main.cpp","$Build/Release/uvlas_aero.lib",
        "$Build/Release/uvlas_core.lib","$prefix/lib/Catch2.lib",'psapi.lib',
        "$prefix/lib/boost_container-vc145-mt-x64-1_92.lib",
        "$prefix/lib/boost_random-vc145-mt-x64-1_92.lib",'/link',"/OUT:$dir/profile.exe","/PDB:$dir/link.pdb")
    & cl.exe @arguments *> "$dir/build.log"
    if($LASTEXITCODE -ne 0){throw "Build failed: $variant"}
}
# Run through tools/dev/msvc.ps1. Pause competing jobs before entering this loop.
# The saved harness fixes its own affinity to the lowest permitted logical CPU.
foreach($round in 0..2){
    foreach($variant in @('before','after')){
        & "$destination/$variant/profile.exe" > "$destination/$variant-$round.txt"
        if($LASTEXITCODE -ne 0){throw "Evaluation failed: $variant/$round"}
    }
}