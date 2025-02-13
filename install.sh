# install R packages
R -e "install.packages(c('BiocManager', 'data.table', 'statmod', 'R.utils'), repos='https://cloud.r-project.org/')"
R -e "BiocManager::install('edgeR')"

# install perl packages
cpanm XML::Compile::Transport::SOAPHTTP XML::Compile::WSDL11 XML::Compile::SOAP11 XML::Compile XML::LibXML::Simple XML::LibXML Log::Log4perl Math::CDF CGI File::Which Data::Dumper File::Copy File::Spec::Functions HTML::PullParser HTML::Template HTML::TreeBuilder JSON Pod::Usage XML::Simple XML::Parser::Expat

# install jbrowse via nodejs
npm install -g @jbrowse/cli
