# Keep Windows CI on the STL whose allocator behavior the input parser models.
$ErrorActionPreference = 'Stop'
if (-not $env:RUNNER_TEMP -or -not $env:GITHUB_ENV) {
    throw 'RUNNER_TEMP and GITHUB_ENV are required.'
}

$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
$installations = @(& $vswhere -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath)
if ($LASTEXITCODE -ne 0) { throw 'Could not enumerate Visual Studio installations.' }

$installation = $installations | Where-Object {
    $versionFile = Join-Path $_ 'VC\Auxiliary\Build\Microsoft.VCToolsVersion.default.txt'
    (Test-Path -LiteralPath $versionFile) -and (Get-Content -LiteralPath $versionFile -Raw).Trim().StartsWith('14.51.')
} | Select-Object -First 1

if (-not $installation) {
    if ($env:GITHUB_ACTIONS -cne 'true') {
        throw 'Install MSVC 14.51 before running this check outside GitHub Actions.'
    }
    # Fixed Visual Studio Build Tools 18.11.3 bootstrapper from Microsoft's release history.
    $hash = 'ea54b67bfac576b6dd097f84eeaf6893780d2a613f90448ff038e8a2777a0afb'
    $url = "https://download.visualstudio.microsoft.com/download/pr/5f474104-723f-4007-bf4c-92e42a23dcc6/$hash/vs_BuildTools.exe"
    $bootstrapper = Join-Path $env:RUNNER_TEMP 'uvlas-vs-BuildTools-18.11.3.exe'
    $installation = Join-Path $env:RUNNER_TEMP 'uvlas-msvc-14.51'
    Invoke-WebRequest -Uri $url -OutFile $bootstrapper -TimeoutSec 120
    if ((Get-FileHash -LiteralPath $bootstrapper -Algorithm SHA256).Hash -ine $hash) {
        throw 'Visual Studio Build Tools bootstrapper hash differs from the pinned release.'
    }
    $options = @('--quiet', '--wait', '--norestart', '--installPath', ('"{0}"' -f $installation),
                 '--add', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64')
    $installer = Start-Process -FilePath $bootstrapper -ArgumentList $options -WindowStyle Hidden -Wait -PassThru
    if ($installer.ExitCode -notin @(0, 3010)) {
        throw "Visual Studio Build Tools install failed with exit code $($installer.ExitCode)."
    }
}

$compilerCheck = Join-Path $env:RUNNER_TEMP 'uvlas-ci-msvc-stl.cpp'
@'
#include <version>

static_assert(_MSC_VER == 1951, "Windows CI requires MSVC 14.51.");
static_assert(_MSVC_STL_VERSION == 145 && _MSVC_STL_UPDATE == 202604L,
              "Windows CI requires the STL revision qualified for input parsing.");
'@ | Set-Content -LiteralPath $compilerCheck -Encoding utf8NoBOM

# Use the same instance for the product and vcpkg dependency builds.
$env:UVLAS_MSVC_INSTALLATION = $installation
try {
    & (Join-Path $PSScriptRoot 'msvc.ps1') cl /nologo /Bv /std:c++23preview /Zs $compilerCheck
    if ($LASTEXITCODE -ne 0) { throw 'The selected MSVC installation failed the toolchain check.' }
} finally {
    Remove-Item -LiteralPath $compilerCheck -Force
}
"UVLAS_MSVC_INSTALLATION=$installation" >> $env:GITHUB_ENV
"VCPKG_VISUAL_STUDIO_PATH=$installation" >> $env:GITHUB_ENV
