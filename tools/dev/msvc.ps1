# Runs one command inside the newest Visual Studio C++ developer environment.
#   pwsh tools/dev/msvc.ps1 cmake --preset windows-msvc
# The Windows presets use Ninja, so the compiler comes from the environment
# rather than from a generator tied to one Visual Studio release.
$ErrorActionPreference = "Stop"
if (-not $env:VSCMD_VER) {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path $vswhere)) { throw "vswhere.exe not found: install Visual Studio with the C++ workload." }
    $install = & $vswhere -latest -prerelease -products * `
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $install) { throw "No Visual Studio installation with the x64 C++ tools was found." }
    $vcvars = Join-Path $install "VC\Auxiliary\Build\vcvars64.bat"
    # The Visual Studio environment points VCPKG_ROOT at its bundled vcpkg.
    # The project uses the vcpkg the caller configured, so keep that value.
    $vcpkgRoot = $env:VCPKG_ROOT
    foreach ($line in (cmd /d /c "`"$vcvars`" >nul 2>&1 && set")) {
        $name, $value = $line -split "=", 2
        if ($name -and $null -ne $value) { Set-Item -Path "Env:$name" -Value $value }
    }
    if (-not $env:VSCMD_VER) { throw "vcvars64.bat did not set up the environment: $vcvars" }
    if ($vcpkgRoot) { $env:VCPKG_ROOT = $vcpkgRoot }
}
if ($args.Count -eq 0) { throw "Usage: msvc.ps1 <command> [arguments...]" }
$command, $rest = $args
& $command @rest
exit $LASTEXITCODE
