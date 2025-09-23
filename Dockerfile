FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

USER root
RUN apt-get update && apt-get install -y tzdata curl gpg

# nodejs packages
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_21.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list

# R packages
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/' | tee -a /etc/apt/sources.list && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# update and install
RUN apt-get update && \
    apt-get install -y libnlopt-dev autoconf automake build-essential ca-certificates cmake cpanminus curl curl fort77 gcc gfortran ghostscript git gnupg gnupg2 gzip libblas-dev libbz2-dev libcurl4-openssl-dev libexpat1-dev libgd-dev libgmp3-dev libgs-dev libgsl-dev libhtml-template-compiled-perl liblapack-dev liblzma-dev libmpc-dev libmpfr-dev libncurses5-dev libopenmpi-dev libssl-dev libtool libxml-libxml-debugging-perl libxml-opml-simplegen-perl libxml2 libxml2-dev libxslt1-dev libxslt1.1 nodejs openmpi-bin python-is-python3 python3 python3-dev python3-pip r-base wget zlib1g-dev && \
    pip3 install deeptools cython --break-system-packages && \
    npm install -g @jbrowse/cli

# R with edgeR
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install('edgeR')" && \
    R -e "install.packages(c('data.table', 'statmod', 'R.utils', 'doParallel', 'foreach', 'nloptr'), repos='https://cloud.r-project.org')"

# Perl modules
RUN cpanm XML::Compile::Transport::SOAPHTTP XML::Compile::WSDL11 XML::Compile::SOAP11 XML::Compile XML::LibXML::Simple XML::LibXML Log::Log4perl Math::CDF CGI File::Which Data::Dumper File::Copy File::Spec::Functions HTML::PullParser HTML::Template HTML::TreeBuilder JSON Pod::Usage XML::Simple XML::Parser::Expat

# MEME suite
RUN wget http://meme-suite.org/meme-software/5.5.6/meme-5.5.6.tar.gz --no-check-certificate && \
    tar zxf meme-5.5.6.tar.gz
WORKDIR /meme-5.5.6
RUN ./configure --prefix=/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt --disable-dependency-tracking && \
    make && \
    make install

# rmats-turbo
RUN apt-get update && \
    apt-get install -y subread bedtools samtools tabix
WORKDIR /
RUN git clone https://github.com/Xinglab/rmats-turbo.git
WORKDIR rmats-turbo
RUN ./build_rmats

# Install the splicekit package
WORKDIR /usr/splicekit
COPY . /usr/splicekit

# Install the package from the current repository
RUN pip install . --break-system-packages

# Install STAR read aligner
WORKDIR /
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz

# environment variables
ENV PATH="$PATH:/rmats-turbo:/meme/bin:/usr/local/bin:/STAR-2.7.11b/bin/Linux_x86_64"

# Run the application
CMD ["splicekit"]
