# SPDX-License-Identifier: Apache-2.0
# Record the outputs of a successful linker invocation before any later assembly.
file(SHA256 "${REVISION_HEADER}" revision_digest)
file(SHA256 "${ARTIFACT}" artifact_digest)
file(WRITE "${OUTPUT}" "${revision_digest}\n${artifact_digest}\n")
