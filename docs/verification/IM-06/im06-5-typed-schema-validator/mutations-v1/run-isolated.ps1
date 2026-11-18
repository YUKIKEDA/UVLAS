param(
    [Parameter(Mandatory=$true)][string]$SourceCommit,
    [string]$Root=(Get-Location).Path,
    [string]$Build='out/issue126-build',
    [string]$CompileDatabase='out/issue126-tidy/compile_commands.json',
    [string]$Output='out/issue126-isolated-mutations-v1'
)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $Root
$sourceHead=& git rev-parse HEAD
if($LASTEXITCODE -ne 0 -or $sourceHead.Trim() -ne $SourceCommit){throw 'Use the fixed source head.'}
& git diff --quiet $SourceCommit -- src/io include/uvlas/io tests/io
if($LASTEXITCODE -ne 0){throw 'Product and test sources must match the fixed commit.'}
$destination=[IO.Path]::GetFullPath((Join-Path $Root $Output))
$allowed=[IO.Path]::GetFullPath((Join-Path $Root 'out'))+[IO.Path]::DirectorySeparatorChar
if(-not $destination.StartsWith($allowed,[StringComparison]::OrdinalIgnoreCase)){throw 'Output escapes workspace out.'}
if(Test-Path -LiteralPath $destination){throw 'Output must be new.'}
New-Item -ItemType Directory -Path "$destination/supply" | Out-Null
$database=Get-Content $CompileDatabase -Raw | ConvertFrom-Json
$command=@($database | Where-Object {$_.file -like '*case_schema.cpp'})[0].command
$found=[regex]::Match($command,'-external:I(?<prefix>[^\s"]+)')
if(-not $found.Success){throw 'No fixed dependency include path.'}
$prefix=Split-Path -Parent $found.Groups['prefix'].Value
$objects=@(Get-ChildItem "$Build/tests/CMakeFiles/uvlas_test_io_case_schema.dir/Debug/io/*.obj")
if($objects.Count -ne 5){throw 'Expected all five fixed normal test objects.'}
foreach($object in $objects){Copy-Item -LiteralPath $object.FullName -Destination "$destination/supply/$($object.Name)"}
foreach($library in @('uvlas_io','uvlas_core')){Copy-Item -LiteralPath "$Build/Debug/$library.lib" -Destination "$destination/supply/$library.lib"}
$cases=@(
    @{Name='control';File='';Filter='[IO4-01]';Edits=@()},
    @{Name='unknown-key-guard';File='case_schema.cpp';Filter='Every structural error*';Edits=@(
        @{Before='context.reject(SchemaReason::UnknownField);';After='member = entry.next; continue;'},
        @{Before='return candidate.name == key;' + "`n" + '                    });';
          After='return candidate.name == key;' + "`n" + '                    });' + "`n" + '                if (field == selected.fields.end()) continue;'}
    )},
    @{Name='exact-tail-16';File='case_schema_exact.cpp';Filter='Exact decimal intervals retain*';Edits=@(
        @{Before='i < std::max(left.size(), right.size())';After='i < std::min(std::max(left.size(), right.size()), std::size_t{16})'}
    )},
    @{Name='shape-plan-omitted';File='case_schema.cpp';Filter='All five frozen root designs*';Edits=@(
        @{Before='shape_pass(context, stack, true);';After='if (context.whole_case && context.case_registry) shape_pass(context, stack, true);'}
    )},
    @{Name='reference-phase-omitted';File='case_schema_exact.cpp';Filter='All five frozen root designs*';Edits=@(
        @{Before='resolve_references(check);';After='if (context.whole_case && context.case_registry) resolve_references(check);'}
    )},
    @{Name='consumer-reference-domain';File='case_schema_exact.cpp';Filter='All expression consumers*';Edits=@(
        @{Before='check.require(domain_kind == "reference", SchemaReason::InvalidExpression);';After='check.require(!domain_kind.empty(), SchemaReason::InvalidExpression);'}
    )},
    @{Name='consumer-time-domain';File='case_schema_exact.cpp';Filter='All expression consumers*';Edits=@(
        @{Before='check.require(domain_kind == "time", SchemaReason::InvalidExpression);';After='check.require(!domain_kind.empty(), SchemaReason::InvalidExpression);'}
    )},
    @{Name='consumer-current-domain';File='case_schema_exact.cpp';Filter='All expression consumers*';Edits=@(
        @{Before='check.require(domain_kind == "current_position_time",';After='check.require(!domain_kind.empty(),'}
    )},
    @{Name='static-window-omitted';File='case_schema_exact.cpp';Filter='Static observations reject*';Edits=@(
        @{Before='if (check.tag(analysis) == "static") {';After='if (check.tag(analysis) == "disabled_static_window") {'}
    )}
)
$summary=@()
foreach($case in $cases){
    $name=$case.Name
    $folder=Join-Path $destination $name
    New-Item -ItemType Directory -Path $folder | Out-Null
    Copy-Item src/io/*.hpp -Destination $folder
    Copy-Item -LiteralPath src/io/case_schema.cpp,src/io/case_schema_exact.cpp -Destination $folder
    $sourceHashes=@{}
    foreach($unit in @('case_schema.cpp','case_schema_exact.cpp')){
        $sourceHashes[$unit]=(Get-FileHash "$folder/$unit" -Algorithm SHA256).Hash.ToLowerInvariant()
    }
    if($case.File){
        $path=Join-Path $folder $case.File
        $original=[IO.File]::ReadAllText($path)
        $crlf=$original.Contains("`r`n")
        $text=$original.Replace("`r`n","`n")
        foreach($edit in $case.Edits){
            $before=$edit.Before;$after=$edit.After
            if(([regex]::Matches($text,[regex]::Escape($before))).Count -ne 1){throw "Unique mutation anchor required: $name"}
            $text=$text.Replace($before,$after)
        }
        if($crlf){$text=$text.Replace("`n","`r`n")}
        [IO.File]::WriteAllText($path,$text)
    }
    $newObjects=@()
    foreach($unit in @('case_schema.cpp','case_schema_exact.cpp')){
        $obj=Join-Path $folder "$unit.obj"
        if($name -ne 'control' -and $unit -ne $case.File){
            Copy-Item -LiteralPath "$destination/control/$unit.obj" -Destination $obj
        }else{
            $arguments=@('/nologo','/TP','/DNOMINMAX','/DUVLAS_NLOHMANN_JSON_ALLOCATION_REPAIR=1',
                '/DWIN32_LEAN_AND_MEAN','/DWIN32','/D_WINDOWS','/EHsc','/Ob0','/Od','/RTC1',
                '/std:c++23preview','/MDd','/Zi','/W4','/permissive-','/utf-8','/Zc:__cplusplus',
                '/fp:strict','/WX',"/I$folder","/I$Root/include","/external:I$prefix/include",'/external:W0',
                "/Fo$obj","/Fd$folder/compile.pdb",'/FS','/c',"$folder/$unit")
            & cl.exe @arguments *> "$folder/$unit.compile.log"
            if($LASTEXITCODE -ne 0){throw "Compile failed: $name/$unit"}
        }
        $newObjects+=$obj
    }
    $exe=Join-Path $folder 'schema-mutation.exe'
    $fixedObjects=@($objects | ForEach-Object {"$destination/supply/$($_.Name)"})
    $link=@('/nologo')+$newObjects+$fixedObjects+@("$destination/supply/uvlas_io.lib",
        "$destination/supply/uvlas_core.lib","$prefix/debug/lib/manual-link/Catch2Maind.lib",
        "$prefix/debug/lib/Catch2d.lib",'psapi.lib','/link',"/OUT:$exe",
        "/PDB:$folder/link.pdb",'/DEBUG','/INCREMENTAL:NO','/SUBSYSTEM:CONSOLE')
    & cl.exe @link *> "$folder/link.log"
    if($LASTEXITCODE -ne 0){throw "Link failed: $name"}
    & $exe $case.Filter --abort *> "$folder/test.log"
    $tested=$LASTEXITCODE
    $log=[IO.File]::ReadAllText("$folder/test.log")
    if($name -eq 'control'){
        if($tested -ne 0 -or $log -notmatch 'All tests passed'){throw 'Control must pass.'}
    }elseif($tested -notin @(1,42) -or $log -notmatch 'FAILED'){
        throw "Expected a Catch assertion failure, not a crash: $name ($tested)"
    }
    $mutatedHashes=@{}
    foreach($unit in @('case_schema.cpp','case_schema_exact.cpp')){
        $mutatedHashes[$unit]=(Get-FileHash "$folder/$unit" -Algorithm SHA256).Hash.ToLowerInvariant()
    }
    $summary+=@{name=$name;source_commit=$SourceCommit;compile_exit=0;link_exit=0;test_exit=$tested;
        filter=$case.Filter;original_hashes=$sourceHashes;mutated_hashes=$mutatedHashes;edits=$case.Edits}
    $summary | ConvertTo-Json -Depth 10 | Set-Content "$destination/summary.json"
    Write-Output "$name compile=0 link=0 test=$tested"
}
