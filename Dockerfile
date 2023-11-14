FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

USER root
RUN yes | unminimize
RUN apt-get update
