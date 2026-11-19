param([Parameter(Mandatory)][string]$Root,[Parameter(Mandatory)][string]$Output,[Parameter(Mandatory)][string]$SourceCommit,[switch]$Resume)
$ErrorActionPreference='Stop'
$Root=(Resolve-Path -LiteralPath $Root).Path
if((Test-Path -LiteralPath $Output)-and !$Resume){throw 'Refusing to overwrite output'}
$FrozenCopy=Join-Path (Split-Path $Output -Parent) 'issue148-fixed-source'
if(!(Test-Path $FrozenCopy)){throw 'Fixed source copy missing'}
$frozen=@{}
foreach($path in (git -C $Root ls-tree -r --name-only $SourceCommit -- src/io src/structure include tests/io)){
 $expected=(git -C $Root rev-parse ($SourceCommit+':'+$path)).Trim()
 $actual=(git -C $Root hash-object -- (Join-Path $FrozenCopy $path)).Trim()
 if($expected-ne $actual){throw "Fixed source mismatch $path"};$frozen[$path]=$actual
}
New-Item -ItemType Directory -Path $Output -Force | Out-Null
$Output=(Resolve-Path -LiteralPath $Output).Path
$build=Join-Path $Root 'out/issue148-build'
$entries=Get-Content (Join-Path $build 'compile_commands.json') -Raw|ConvertFrom-Json
$files=@('si_case.cpp','si_case_storage.cpp','si_case_catalog.cpp','si_quantity.cpp','si_exact.cpp','si_exact_center.cpp')
$supply=Join-Path $Output 'supply';New-Item -ItemType Directory -Path $supply -Force | Out-Null
foreach($lib in @('uvlas_io','uvlas_structure','uvlas_core')){Copy-Item (Join-Path $PSScriptRoot "issue148-mutations/supply/$lib.lib") $supply}
$normalObjects=@()
foreach($test in @('units','ownership','fragments','fields')){$name="si_case_${test}_test.cpp.obj";Copy-Item (Join-Path $PSScriptRoot "issue148-mutations/supply/$name") $supply;$normalObjects+=(Join-Path $supply $name)}
$prefix='C:/Users/yuuuu/.codex/worktrees/uvlas-im00-9/UVLAS/out/issue99-vcpkg-installed/x64-windows'
$cases=@(
 @('control','','','','normal'),
 @('length-factor-inverse','si_exact.cpp','si::divide(si::integer(length_n[l]), si::integer(length_d[l]), c.work)','si::divide(si::integer(length_d[l]), si::integer(length_n[l]), c.work)','normal'),
 @('degree-pi-factor-replaced','si_exact.cpp','si::programmed_power(pi, recipe.pi_power, c.program.values()[0], c.work)','si::integer(1)','normal'),
 @('two-word-low-dropped','si_exact_center.cpp','ScaledTwoWord{h.word, l.word, 0, q}','ScaledTwoWord{h.word, 0.0, 0, q}','normal'),
 @('foreign-catalog-checks-disabled','si_case_catalog.cpp','if (!detail::SiFactory::owns_key(*source, row.key))','if (false && !detail::SiFactory::owns_key(*source, row.key))','normal'),
 @('work-meter-disabled','si_interval.hpp','used += amount;','used += amount * 0;','interval'),
 @('division-directed-rounding','si_interval.hpp','if (enlarged % denominator != 0 && upper != negative)','if (enlarged % denominator != 0 && upper == negative)','interval')
)
$results=@();$controlObjects=@{}
if($Resume){$old=Get-Content (Join-Path $Output 'summary.json') -Raw|ConvertFrom-Json -AsHashtable;if($old.source_commit-ne $SourceCommit){throw 'Resume source mismatch'};$results=@($old.cases);foreach($src in $files){$controlObjects[$src]=Join-Path $Output "control/$src.obj"}}
function Compile-Source($source,$directory,$test){
 $entry=@($entries|Where-Object {(Split-Path $_.file -Leaf)-eq $source -and $_.output -like '*Debug*' -and $(if($test){$_.output -like '*uvlas_test_io_si_interval.dir*'}else{$_.output -like '*uvlas_io.dir*'})})[0]
 if(!$entry){throw "Compile entry missing $source"}
 $tokens=@($entry.command -split '\s+');$compiler=$tokens[0]
 $args=@($tokens|Select-Object -Skip 1|Where-Object {$_ -notlike '/Fo*' -and $_ -notlike '/Fd*' -and $_ -notlike '-DCMAKE_INTDIR*' -and $_ -ne $entry.file}|ForEach-Object {$_.Replace($Root,$directory).Replace($Root.Replace('\','/'),$directory.Replace('\','/'))})
 $object=Join-Path $directory "$source.obj";$relative=if($test){'tests/io/'}else{'src/io/'}
 $args+=@("/Fo$object","/Fd$directory/compile.pdb",(Join-Path $directory ($relative+$source)))
 & $compiler @args *> (Join-Path $directory "$source.compile.txt")
 if($LASTEXITCODE-ne 0){throw "Compile failed $source in $directory"};return $object
}
function Link-Run($objects,$directory,$name){
 $exe=Join-Path $directory "$name.exe"
 $args=@('/nologo')+$objects+@((Join-Path $supply 'uvlas_io.lib'),(Join-Path $supply 'uvlas_structure.lib'),(Join-Path $supply 'uvlas_core.lib'),"$prefix/debug/lib/manual-link/Catch2Maind.lib","$prefix/debug/lib/Catch2d.lib",'psapi.lib','/link',"/OUT:$exe","/PDB:$directory/$name.pdb",'/DEBUG','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE')
 & cl.exe @args *> (Join-Path $directory "$name.link.txt")
 if($LASTEXITCODE-ne 0){throw "Link failed $directory $name"}
 & $exe --abort --reporter compact *> (Join-Path $directory "$name.test.txt")
 return [int]$LASTEXITCODE
}
foreach($case in $cases){
 $name,$file,$before,$after,$kind=$case;if($results.name -contains $name){continue};$d=Join-Path $Output $name
 New-Item -ItemType Directory -Path $d | Out-Null
 Copy-Item (Join-Path $FrozenCopy 'include') $d -Recurse
 New-Item -ItemType Directory -Path (Join-Path $d 'src'),(Join-Path $d 'tests') | Out-Null
 Copy-Item (Join-Path $FrozenCopy 'src/io') (Join-Path $d 'src') -Recurse
 Copy-Item (Join-Path $FrozenCopy 'src/structure') (Join-Path $d 'src') -Recurse
 Copy-Item (Join-Path $FrozenCopy 'tests/io') (Join-Path $d 'tests') -Recurse
 if($file){$path=Join-Path $d "src/io/$file";$s=[IO.File]::ReadAllText($path);if(([regex]::Matches($s,[regex]::Escape($before))).Count-ne 1){throw "Anchor not unique $name"};[IO.File]::WriteAllText($path,$s.Replace($before,$after),[Text.UTF8Encoding]::new($false))}
 if($name-eq 'foreign-catalog-checks-disabled'){$path=Join-Path $d 'src/io/si_case_catalog.cpp';$s=[IO.File]::ReadAllText($path);$a='if (!(row.key == entry.key))';if(([regex]::Matches($s,[regex]::Escape($a))).Count-ne 1){throw 'Second generation guard anchor'};[IO.File]::WriteAllText($path,$s.Replace($a,'if (false && !(row.key == entry.key))'),[Text.UTF8Encoding]::new($false))}
 if($kind-eq 'interval'){$objs=@((Compile-Source 'si_case_interval_test.cpp' $d $true));$exit=Link-Run $objs $d 'interval'}else{
  $objs=@();foreach($src in $files){if($name-ne 'control' -and $src-ne $file){$objs+=$controlObjects[$src]}else{$o=Compile-Source $src $d $false;$objs+=$o;if($name-eq 'control'){$controlObjects[$src]=$o}}}
  $exit=Link-Run ($objs+$normalObjects) $d 'normal'
  if($name-eq 'control'){$interval=@((Compile-Source 'si_case_interval_test.cpp' $d $true));$i=Link-Run $interval $d 'interval';if($i-ne 0){throw 'Interval control failed'}}
 }
 if(($name-eq 'control' -and $exit-ne 0)-or($name-ne 'control' -and $exit-notin @(1,42))){throw "Unexpected Catch $name $exit"}
 $results+=@{name=$name;target=$kind;compile_exit=0;link_exit=0;test_exit=$exit}
 [IO.File]::WriteAllText((Join-Path $Output 'summary.json'),((@{source_commit=$SourceCommit;source_git_blobs=$frozen;cases=$results}|ConvertTo-Json -Depth 8).Replace("`r`n","`n")+"`n"),[Text.UTF8Encoding]::new($false))
 Write-Output "$name compile0 link0 Catch$exit"
}
$global:LASTEXITCODE=0
exit 0