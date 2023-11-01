FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

USER root
RUN yes | unminimize
RUN echo '8'| apt-get install -y tzdata curl gpg

# nodejs packages
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_21.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list

# R packages
#RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' | tee -a /etc/apt/sources.list
#RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# update and install
RUN apt-get update
RUN apt-get install -y autoconf automake build-essential ca-certificates cmake cpanminus curl curl fort77 gcc gfortran ghostscript git gnupg gnupg2 gzip libblas-dev libbz2-dev libcurl4-openssl-dev libexpat1-dev libgd-dev libgmp3-dev libgs-dev libgsl-dev libhtml-template-compiled-perl liblapack-dev liblzma-dev libmpc-dev libmpfr-dev libncurses5-dev libopenmpi-dev libssl-dev libtool libxml-libxml-debugging-perl libxml-opml-simplegen-perl libxml2 libxml2-dev libxslt1-dev libxslt1.1 nodejs openmpi-bin python-is-python3 python3 python3-dev python3-pip r-base wget zlib1g-dev
RUN pip3 install deeptools cython==0.29.36
RUN npm install -g @jbrowse/cli

# R with edgeR
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('statmod')"

# Perl modules
RUN cpanm XML::Compile::Transport::SOAPHTTP XML::Compile::WSDL11 XML::Compile::SOAP11 XML::Compile XML::LibXML::Simple XML::LibXML Log::Log4perl Math::CDF CGI File::Which Data::Dumper File::Copy File::Spec::Functions HTML::PullParser HTML::Template HTML::TreeBuilder JSON Pod::Usage XML::Simple XML::Parser::Expat

# MEME suite
RUN wget http://meme-suite.org/meme-software/5.5.4/meme-5.5.4.tar.gz
RUN tar zxf meme-5.5.4.tar.gz
WORKDIR /meme-5.5.4
RUN ./configure --prefix=/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt --disable-dependency-tracking
RUN make
RUN make install

# rmats-turbo
RUN apt-get install -y subread bedtools samtools tabix
WORKDIR /
RUN git clone https://github.com/Xinglab/rmats-turbo.git
WORKDIR rmats-turbo
RUN ./build_rmats

# environment variables
ENV PATH="$PATH:/rmats-turbo:/meme/bin:/usr/local/bin"
