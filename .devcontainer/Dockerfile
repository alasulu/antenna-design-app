# Pinned to the digest last built and checked (codespace workflow): a floating tag
# changed under us once already. Bump it deliberately, with that workflow's run.
FROM mcr.microsoft.com/devcontainers/python:1-3.12-bookworm@sha256:7876580dc67fd460fd962f004cbeb480027e9bbc0657096f1087db11f9eaff39
# The image lists Yarn's apt source, whose signing key no longer verifies; every
# apt-get update then fails, the virtual desktop's install included. Nothing here
# needs Yarn.
RUN rm -f /etc/apt/sources.list.d/yarn.list
