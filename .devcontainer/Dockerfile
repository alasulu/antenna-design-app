FROM mcr.microsoft.com/devcontainers/python:1-3.12-bookworm
# The image lists Yarn's apt source, whose signing key no longer verifies; every
# apt-get update then fails, the virtual desktop's install included. Nothing here
# needs Yarn.
RUN rm -f /etc/apt/sources.list.d/yarn.list
