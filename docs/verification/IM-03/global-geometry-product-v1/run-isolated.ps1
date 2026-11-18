param([Parameter(Mandatory)][string]$Root,[Parameter(Mandatory)][string]$Output,[Parameter(Mandatory)][string]$SourceCommit)
$ErrorActionPreference='Stop'
$Root=(Resolve-Path -LiteralPath $Root).Path
if($Root.Contains(' ')){throw 'This fixed MSVC evidence producer requires a space-free checkout path'}
if((git -C $Root ls-files --others --exclude-standard src/aero include)){throw 'Untracked source would contaminate the fixed copy'}
$frozen=@{}
foreach($path in (git -C $Root ls-files src/aero include tests/aero/global_geometry_test.cpp tests/aero/global_geometry_fixture.hpp tests/aero/global_kinematics_fixture.hpp)){
    $expected=(git -C $Root rev-parse ($SourceCommit+':'+$path)).Trim()
    if($LASTEXITCODE-ne 0){throw "Missing source at commit: $path"}
    $actual=(git -C $Root hash-object -- (Join-Path $Root $path)).Trim()
    if($actual-ne $expected){throw "Source differs from fixed commit: $path"}
    $frozen[$path]=$expected
}
if(Test-Path -LiteralPath $Output){throw 'Refusing to overwrite evidence output'}
New-Item -ItemType Directory -Path $Output | Out-Null
$Output=(Resolve-Path -LiteralPath $Output).Path
$build=Join-Path $Root 'out/issue147-build'
$entries=Get-Content -LiteralPath (Join-Path $Root 'out/issue147-tidy-db/compile_commands.json') -Raw|ConvertFrom-Json
$files=@('global_geometry.cpp','global_bvh.cpp','global_intersections.cpp','alias_kinematics.cpp')
$supply=Join-Path $Output 'supply'
New-Item -ItemType Directory -Path $supply | Out-Null
Copy-Item -LiteralPath (Join-Path $build 'Debug/uvlas_aero.lib'),(Join-Path $build 'Debug/uvlas_core.lib'),(Join-Path $build 'tests/CMakeFiles/uvlas_test_aero_global_geometry.dir/Debug/aero/global_geometry_test.cpp.obj') -Destination $supply
$prefix=''
foreach($token in ($entries[0].command -split '\s+')){if($token.StartsWith('-external:I')){$prefix=Split-Path -Parent $token.Substring(11)}}
if(-not $prefix){throw 'External dependency prefix missing'}
$cases=@(
    @('control','','',''),
    @('second-pass','global_geometry.cpp','const auto traversal_budget = times(2, times(v, node_cost));','const auto traversal_budget = times(1, times(v, node_cost));'),
    @('alias-basis','global_geometry.cpp','raw_incidence * basis != facet.incidence[k]','raw_incidence != facet.incidence[k]'),
    @('TE-sign','global_geometry.cpp','static_cast<std::int8_t>(sign * incidence)','static_cast<std::int8_t>(-sign * incidence)'),
    @('internal-TE','global_geometry.cpp','if (count != 1)','if (count == 0)'),
    @('shared-point-skip','global_intersections.cpp','const auto& fb = g.facets.values()[static_cast<std::size_t>(ib)];','const auto& fb = g.facets.values()[static_cast<std::size_t>(ib)]; for(const auto x:fa.points)for(const auto y:fb.points)if(x==y)return false;'),
    @('open-box','global_bvh.cpp','if (compare(al, bh, r) > 0)','if (compare(al, bh, r) >= 0)'),
    @('full-J','alias_kinematics.cpp','!equal_scalar(a.value, b.value, run)','(equal_scalar(a.value, b.value, run) && false)'),
    @('explicit-time','alias_kinematics.cpp','!equal_scalar(p.time[k], owner.time[k], run)','(equal_scalar(p.time[k], owner.time[k], run) && false)'),
    @('velocity-assertion','alias_kinematics.cpp','if (point.asserted)','if (point.asserted && point.id == UINT64_MAX)'),
    @('epoch','global_geometry.cpp','return MaterialReleaseKey{*canonical, epoch};','return MaterialReleaseKey{*canonical, ReleaseEpoch{epoch.value()*0}};')
)
$results=@()
$controlObjects=@{}
foreach($case in $cases){
    $name,$file,$before,$after=$case
    $directory=Join-Path $Output $name
    New-Item -ItemType Directory -Path (Join-Path $directory 'src/aero') -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $Root 'include') -Destination $directory -Recurse
    Get-ChildItem -LiteralPath (Join-Path $Root 'src/aero') -File | Copy-Item -Destination (Join-Path $directory 'src/aero')
    if($file){$path=Join-Path $directory "src/aero/$file";$text=[IO.File]::ReadAllText($path);if(([regex]::Matches($text,[regex]::Escape($before))).Count-ne 1){throw "Mutation anchor count: $name"};[IO.File]::WriteAllText($path,$text.Replace($before,$after),[Text.UTF8Encoding]::new($false))}
    $objects=@()
    foreach($source in $files){
        if($name-ne 'control' -and $source-ne $file){$objects+=$controlObjects[$source];continue}
        $entry=@($entries|Where-Object {(Split-Path $_.file -Leaf)-eq $source})[0]
        $tokens=@($entry.command -split '\s+')
        $compiler=$tokens[0]
        $arguments=@($tokens|Select-Object -Skip 1|Where-Object {$_ -notlike '/Fo*' -and $_ -notlike '/Fd*' -and $_ -notlike '-DCMAKE_INTDIR*' -and $_ -ne $entry.file}|ForEach-Object {$_.Replace($Root,$directory)})
        $object=Join-Path $directory "$source.obj"
        $arguments+=@("/Fo$object","/Fd$directory/compile.pdb",(Join-Path $directory "src/aero/$source"))
        & $compiler @arguments *> (Join-Path $directory "$source.compile.txt")
        if($LASTEXITCODE-ne 0){throw "Compile failed: $name / $source"}
        $objects+=$object
        if($name-eq 'control'){$controlObjects[$source]=$object}
    }
    $exe=Join-Path $directory 'test.exe'
    $link=@('/nologo')+$objects+@((Join-Path $supply 'global_geometry_test.cpp.obj'),(Join-Path $supply 'uvlas_aero.lib'),(Join-Path $supply 'uvlas_core.lib'),"$prefix/debug/lib/manual-link/Catch2Maind.lib","$prefix/debug/lib/Catch2d.lib",'psapi.lib','/link',"/OUT:$exe","/PDB:$directory/link.pdb",'/DEBUG','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE')
    & cl.exe @link *> (Join-Path $directory 'link.txt')
    if($LASTEXITCODE-ne 0){throw "Link failed: $name"}
    & $exe '[UW-D07]' --abort --reporter compact *> (Join-Path $directory 'test.txt')
    $test=$LASTEXITCODE
    if(($name-eq 'control' -and $test-ne 0)-or($name-ne 'control' -and $test-notin @(1,42))){throw "Unexpected Catch exit $name : $test"}
    $results+=@{name=$name;compile_exit=0;link_exit=0;test_exit=$test;filter='[UW-D07]'}
    $json=(@{source_commit=$SourceCommit;source_git_blobs=$frozen;cases=$results}|ConvertTo-Json -Depth 5).Replace("`r`n","`n")+"`n"
    [IO.File]::WriteAllText((Join-Path $Output 'summary.json'),$json,[Text.UTF8Encoding]::new($false))
    Write-Output "$name compile=0 link=0 Catch=$test"
}
