FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

USER root
RUN yes | unminimize
RUN apt-get update
RUN echo '8'| apt-get install -y tzdata
RUN apt-get install -y wget gcc python3.8 curl git libgsl-dev cmake gfortran fort77 libblas-dev liblapack-dev python-is-python3 gnupg2 build-essential gzip python3-pip python3-dev zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-openssl-dev libssl-dev libopenmpi-dev openmpi-bin ghostscript libgs-dev libgd-dev libexpat1-dev zlib1g-dev libxml2-dev autoconf automake libtool libhtml-template-compiled-perl libxml-opml-simplegen-perl libxml-libxml-debugging-perl
RUN pip3 install deeptools cython==0.29.36
RUN curl -fsSL https://deb.nodesource.com/setup_14.x | bash -
RUN apt-get install -y nodejs
RUN npm install -g @jbrowse/cli

RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' | tee -a /etc/apt/sources.list
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# R with edgeR
RUN apt-get update
RUN apt-get install -y r-base
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('statmod')"

RUN apt-get install -y libxml2 libxml2-dev libxslt1.1 libxslt1-dev libcurl4-gnutls-dev libexpat1-dev zlib1g-dev libgmp3-dev libmpfr-dev libmpc-dev cpanminus
RUN cpanm XML::Compile::Transport::SOAPHTTP XML::Compile::WSDL11 XML::Compile::SOAP11 XML::Compile XML::LibXML::Simple XML::LibXML Log::Log4perl Math::CDF CGI File::Which Data::Dumper File::Copy File::Spec::Functions HTML::PullParser HTML::Template HTML::TreeBuilder JSON Pod::Usage XML::Simple XML::Parser::Expat

# MEME suite
RUN wget http://meme-suite.org/meme-software/5.5.3/meme-5.5.3.tar.gz
RUN tar zxf meme-5.5.3.tar.gz
WORKDIR /meme-5.5.3
RUN ./configure --prefix=/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt --disable-dependency-tracking
RUN make
RUN make install

# rmats-turbo
RUN apt-get install -y subread bedtools samtools tabix
RUN pip3 install pybio
WORKDIR /
RUN git clone https://github.com/Xinglab/rmats-turbo.git
WORKDIR rmats-turbo
RUN ./build_rmats

# environment variables
ENV PATH="$PATH:/rmats-turbo:/meme/bin:/usr/local/bin"
#ENV PYTHONPATH="$PYTHONPATH:/meme-5.5.3/scripts"
#ENV MEME_BIN_DIR="/meme-5.5.3/scripts"
