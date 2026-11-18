param([string]$Root = (Get-Location).Path, [string]$OnlyName = '')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $Root
$base=Join-Path $Root 'out/p1-mutations-final'
New-Item -ItemType Directory -Path $base -Force | Out-Null
$entry=(Get-Content out/issue127-tidy-db/compile_commands.json -Raw | ConvertFrom-Json)[0]
$include=@($entry.arguments | Where-Object {$_ -like '-external:I*'})[0]
$prefix=Split-Path -Parent $include.Substring('-external:I'.Length)
$build=Join-Path $Root 'out/issue127-build'
$cases=@(
    @('control','','','','[UI-D02]'),
    @('low-word','p1_induction.cpp','        out += temporary;','        out += temporary * 0;','P1 huge common*'),
    @('sheet','p1_program_private.hpp','auto velocity = cross(cross(n, gradient), sheet);','auto velocity = cross(cross(n, gradient), sheet) * p.constant(0);','P1 all 126*'),
    @('radius-seed','p1_program_private.hpp','node.left < 13','node.left < 12','P1 all 126*'),
    @('exterior-normal','p1_program_private.hpp','normal = d * jac * 2 / delta;','normal = d * jac * 0 / delta;','P1 all 126*'),
    @('boundary-sign','p1_program_private.hpp','velocity = velocity + line[edge] * strength;','velocity = velocity - line[edge] * strength;','P1 all 126*'),
    @('far-qualification','p1_induction.cpp','if (compare(delta.upper, threshold) > 0)','if (compare(delta.upper, threshold) > 1)','P1 far qualification rejects*'),
    @('half-width','p1_induction.cpp','output_scale) - 2146 + scale','output_scale) - 2148 + scale','P1 exact qualification*'),
    @('self-zero','p1_induction.cpp','p1::scale(evaluated.components[jet], -scale)','p1::scale(std::holds_alternative<VerifiedIncidence>(needed(packet.qualification())) ? p1::integer(0) : evaluated.components[jet], -scale)','P1 all 126*'),
    @('abs-differential','p1_program_private.hpp','out[j] = b[j];','out[j] = b[j].upper.mantissa < 0 ? negate(b[j], work) : b[j];','P1 all 126*'),
    @('scale-differential','p1_program_private.hpp','scale(integer(1), seed_power)','scale(integer(1), seed_power * 0)','P1 huge common*')
)
$summary=if($OnlyName){@(Get-Content (Join-Path $base 'summary.json') -Raw | ConvertFrom-Json)}else{@()}
foreach($case in $cases) {
    $name,$file,$before,$after,$filter=$case
    if($OnlyName -and $OnlyName -ne $name){continue}
    $dir=Join-Path $base $name
    if(Test-Path $dir){throw "Output exists: $dir"}
    New-Item -ItemType Directory -Path $dir | Out-Null
    Copy-Item src/aero/*.hpp -Destination $dir
    Copy-Item src/aero/p1_induction.cpp -Destination $dir
    if($file){
        $path=Join-Path $dir $file
        $text=[IO.File]::ReadAllText($path)
        if(([regex]::Matches($text,[regex]::Escape($before))).Count -ne 1){throw "Replacement count: $name"}
        [IO.File]::WriteAllText($path,$text.Replace($before,$after))
    }
    $obj=Join-Path $dir 'p1.obj'
    $exe=Join-Path $dir 'p1.exe'
    $arguments=@($entry.arguments | Select-Object -Skip 1 | Where-Object {$_ -notlike '/Fo*' -and $_ -notlike '/Fd*' -and $_ -ne $entry.file})
    $arguments+=@("/Fo$obj","/Fd$dir/compile.pdb",(Join-Path $dir 'p1_induction.cpp'))
    & $entry.arguments[0] @arguments *> (Join-Path $dir 'compile.log')
    $compile=$LASTEXITCODE
    if($compile -ne 0){throw "Compile failed $name"}
    $link=@('/nologo',$obj,"$build/tests/CMakeFiles/uvlas_test_aero_p1_induction.dir/Debug/aero/p1_induction_test.cpp.obj","$build/tests/CMakeFiles/uvlas_test_aero_p1_induction.dir/Debug/aero/p1_area_reference_test.cpp.obj","$build/Debug/uvlas_aero.lib","$build/Debug/uvlas_core.lib","$prefix/debug/lib/manual-link/Catch2Maind.lib","$prefix/debug/lib/Catch2d.lib","$prefix/debug/lib/boost_container-vc145-mt-gd-x64-1_92.lib","$prefix/debug/lib/boost_random-vc145-mt-gd-x64-1_92.lib",'psapi.lib','/link',"/OUT:$exe","/PDB:$dir/link.pdb",'/DEBUG','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE')
    & cl.exe @link *> (Join-Path $dir 'link.log')
    $linked=$LASTEXITCODE
    if($linked -ne 0){throw "Link failed $name"}
    & $exe $filter --abort *> (Join-Path $dir 'test.log')
    $tested=$LASTEXITCODE
    if(($name -eq 'control' -and $tested -ne 0) -or ($name -ne 'control' -and $tested -notin @(1,42))){throw "Unexpected test exit $name : $tested"}
    $summary+=@{name=$name; compile_exit=$compile;link_exit=$linked;test_exit=$tested;filter=$filter}
    $summary | ConvertTo-Json -Depth 8 | Set-Content (Join-Path $base 'summary.json')
    Write-Output "$name compile=$compile link=$linked test=$tested"
}
