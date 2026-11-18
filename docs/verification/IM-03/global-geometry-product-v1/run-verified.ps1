param([Parameter(Mandatory)][string]$Root,[Parameter(Mandatory)][string]$Output,
      [Parameter(Mandatory)][string]$SourceCommit,[switch]$ValidateExisting)
$ErrorActionPreference='Stop'
try {
    if(-not $ValidateExisting){& (Join-Path $PSScriptRoot 'run-isolated.ps1') -Root $Root -Output $Output -SourceCommit $SourceCommit}
    $summary=Get-Content -LiteralPath (Join-Path $Output 'summary.json') -Raw|ConvertFrom-Json
    if($summary.source_commit-ne $SourceCommit){throw 'Source identity mismatch'}
    $expected=@('control','second-pass','alias-basis','TE-sign','internal-TE','shared-point-skip','open-box','full-J','explicit-time','velocity-assertion','epoch')
    if($summary.cases.Count-ne $expected.Count){throw 'Incomplete mutation results'}
    $seen=@{}
    foreach($case in $summary.cases){
        if($case.name-notin $expected -or $seen.ContainsKey($case.name)){throw 'Unknown or duplicated mutation'}
        $seen[$case.name]=$true
        foreach($property in @('compile_exit','link_exit','test_exit')){if($case.$property-isnot [long]){throw 'Exit status must be an integer'}}
        if($case.compile_exit-ne 0 -or $case.link_exit-ne 0){throw 'A compile/link failure cannot qualify a mutation'}
        if(($case.name-eq 'control' -and $case.test_exit-ne 0)-or($case.name-ne 'control' -and $case.test_exit-ne 42)){throw 'Unexpected Catch outcome'}
        foreach($log in @('link.txt','test.txt')){if(-not(Test-Path -LiteralPath (Join-Path (Join-Path $Output $case.name) $log))){throw 'Missing execution log'}}
    }
    Write-Output 'Verified: control passed; 10 mutations compiled, linked, and failed Catch assertions.'
    # Expected negative Catch statuses must not become the producer's exit code.
    exit 0
} catch {
    Write-Error $_
    exit 1
}
