# Runs the real wrapper in isolated children without changing the caller's VS environment.
$ErrorActionPreference = 'Stop'
$wrapper = Join-Path $PSScriptRoot 'msvc.ps1'
$scratch = Join-Path ([IO.Path]::GetTempPath()) ([Guid]::NewGuid().ToString('N'))
[void][IO.Directory]::CreateDirectory($scratch)
try {
    $probe = Join-Path $scratch 'argument probe.ps1'
    $probeSource = @'
@{ arguments = @($args) } | ConvertTo-Json -Compress
exit ([int]$env:UVLAS_PROBE_EXIT)
'@
    [IO.File]::WriteAllText($probe, $probeSource, [Text.UTF8Encoding]::new($false))
    $cases = @(
        @{ Arguments = @(); ExitCode = 0 },
        @{ Arguments = @('run.py'); ExitCode = 0 },
        @{ Arguments = @('first', 'second', 'third'); ExitCode = 0 },
        @{ Arguments = @('with spaces', '', 'last'); ExitCode = 0 },
        @{ Arguments = @('run.py'); ExitCode = 23 }
    )
    foreach ($case in $cases) {
        $start = [Diagnostics.ProcessStartInfo]::new((Get-Process -Id $PID).Path)
        $start.UseShellExecute = $false
        $start.RedirectStandardOutput = $true
        $start.RedirectStandardError = $true
        $start.Environment['VSCMD_VER'] = 'isolated-argument-probe'
        $start.Environment['UVLAS_PROBE_EXIT'] = [string]$case.ExitCode
        foreach ($argument in @('-NoProfile', '-File', $wrapper, $probe) + $case.Arguments) {
            $start.ArgumentList.Add($argument)
        }
        $child = [Diagnostics.Process]::Start($start)
        try {
            if (-not $child.WaitForExit(10000)) {
                $child.Kill($true)
                throw 'Wrapper probe timed out.'
            }
            $stdout = $child.StandardOutput.ReadToEnd()
            $stderr = $child.StandardError.ReadToEnd()
            if ($child.ExitCode -ne $case.ExitCode) { throw "Wrong exit code: $($child.ExitCode); $stderr" }
            $actual = @((ConvertFrom-Json $stdout).arguments)
            if ($actual.Count -ne $case.Arguments.Count) { throw "Wrong argument count: $stdout" }
            for ($index = 0; $index -lt $actual.Count; ++$index) {
                if ($actual[$index] -cne $case.Arguments[$index]) { throw "Wrong argument at ${index}: $stdout" }
            }
        }
        finally { $child.Dispose() }
    }
    Write-Output 'Passed 5 MSVC argument/exit-code probes.'
}
finally { Remove-Item -LiteralPath $scratch -Recurse -Force }
