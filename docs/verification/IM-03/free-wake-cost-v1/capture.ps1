param(
    [Parameter(Mandatory)][string] $Executable,
    [Parameter(Mandatory)][string] $Output,
    [int] $Cpu = 0
)

$ErrorActionPreference = 'Stop'
$binary = (Resolve-Path -LiteralPath $Executable).Path
New-Item -ItemType Directory -Path $Output -Force | Out-Null
$run = (Resolve-Path -LiteralPath $Output).Path
$hostProcess = Get-Process -Id $PID
$originalAffinity = $hostProcess.ProcessorAffinity
try {
    $hostProcess.ProcessorAffinity = [IntPtr]([int64]1 -shl $Cpu)
    foreach ($panels in @(2, 4, 8)) {
        $log = Join-Path $run "wake-$panels.log"
        & $binary "[.wake-cost-$panels]" --rng-seed 23 *> $log
        if ($LASTEXITCODE) {
            Get-Content -LiteralPath $log -Tail 24
            exit $LASTEXITCODE
        }
    }
    python -X utf8 (Join-Path $PSScriptRoot 'extract.py') $run
    exit $LASTEXITCODE
} finally {
    $hostProcess.ProcessorAffinity = $originalAffinity
}
