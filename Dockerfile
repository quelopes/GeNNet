FROM centos:7

MAINTAINER Raquel Lopes Costa "quelopes@gmail.com"
EXPOSE 3838 7474 8787

# =============
# --- Linux ---
# =============
RUN yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm && \
    yum -y update && \
    yum -y install libcurl-devel libxml2-devel openssl-devel libpng-devel R wget && \
    yum -y clean all

# ======================================
# --- INSTALL BIOCONDUCTOR AND RNEO4J---
# ======================================
RUN R -e "source(\"https://bioconductor.org/biocLite.R\"); biocLite()" && \
    R -e "source(\"https://bioconductor.org/biocLite.R\"); biocLite(c(\"GOstats\",\"hgu133plus2cdf\",\"hgu133acdf\",\"hgu133a2cdf\",\"hugene10stv1cdf\",\"affy\",\"impute\",\"Biobase\",\"limma\",\"org.Mmu.eg.db\",\"org.Mm.eg.db\",\"org.Rn.eg.db\",\"genefilter\",\"org.Hs.eg.db\",\"ggplot2\",\"igraph\",\"VennDiagram\",\"gplots\", \"fpc\",\"stringr\",\"WGCNA\",\"dynamicTreeCut\",\"frma\"))" && \
    R -e "install.packages('RNeo4j',repos='https://cran.rstudio.com/', clean=TRUE)"

# =======================
# --- INSTALL RSTUDIO ---
# =======================
RUN VER=$(wget --no-check-certificate -qO- https://s3.amazonaws.com/rstudio-server/current.ver) && \
    wget https://download2.rstudio.org/rstudio-server-rhel-${VER}-x86_64.rpm && \
    yum -y install --nogpgcheck rstudio-server-rhel-${VER}-x86_64.rpm && \
    rm rstudio-server-rhel-${VER}-x86_64.rpm && \
    useradd -m rstudio && \
    echo "rstudio:rstudio" | chpasswd && \
    yum -y clean all

# =====================
# --- INSTALL NEO4J ---
# =====================
RUN wget -O neo4j.tar.gz https://neo4j.com/artifact.php?name=neo4j-community-3.0.6-unix.tar.gz && \
    cd /usr/local; tar xvfz /neo4j.tar.gz; ln -s neo4j-community-3.0.6 neo4j && \
    rm /neo4j.tar.gz && \
    sed -i -e "s/#dbms.connector.http.address/dbms.connector.http.address/" /usr/local/neo4j/conf/neo4j.conf


# =================================
# --- BUILD BACKGROUND DATABASE ---
# =================================
ADD import.cql /var/tmp/
ADD DataGraph/data.csv //usr/local/neo4j/data/databases/graph.db/import/
ADD DataGraph/PPI.csv //usr/local/neo4j/data/databases/graph.db/import/
ADD DataGraph/TFs.csv //usr/local/neo4j/data/databases/graph.db/import/
ADD DataGraph/taxon.csv //usr/local/neo4j/data/databases/graph.db/import/

RUN /usr/local/neo4j/bin/neo4j start && \
    sleep 20 && \
    curl -H "Content-Type: application/json" -X POST -d '{"password":"graph"}' -u neo4j:neo4j http://localhost:7474/user/neo4j/password && \
    /usr/local/neo4j/bin/neo4j stop && \
    /usr/local/neo4j/bin/neo4j-shell -v -path /usr/local/neo4j/data/databases/graph.db -config /usr/local/neo4j/conf/neo4j.conf -file /var/tmp/import.cql && \
    rm -rf /usr/local/neo4j/data/databases/graph.db/neostore.transaction.db.*

# =============
# --- Shiny ---
# =============
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/centos-6.3/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/centos-6.3/x86_64/shiny-server-$VERSION-rh6-x86_64.rpm" -O ss-latest.rpm && \
    yum -y install ss-latest.rpm && \
    mkdir /usr/share/doc/R-3.3.1/html/ && \
    rm -f version.txt ss-latest.rpm && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cran.rstudio.com/', clean=TRUE)" && \
    cp -R /usr/lib64/R/library/shiny/examples/* /srv/shiny-server/ && \
    yum -y clean all
RUN R -e "install.packages(c('shinydashboard','shiny','shinythemes','RNeo4j','visNetwork','ggplot2','data.table','networkD3','igraph','shinyBS','RColorBrewer','devtools','d3heatmap'), repos='https://cran.rstudio.com/', clean=TRUE)"

RUN wget -O RDataTracker.tar.gz http://harvardforest.fas.harvard.edu/data/p09/hf091/hf091-01-RDataTracker_2.24.0.tar.gz && \
    R CMD INSTALL RDataTracker.tar.gz && \
    rm RDataTracker.tar.gz

COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod 755 /usr/bin/shiny-server.sh
COPY GennetShiny /srv/shiny-server/gennet


# =========================
# --- BUILD THE MODULES ---
# =========================
USER rstudio
RUN mkdir /home/rstudio/Module-A
RUN mkdir /home/rstudio/Annotation
RUN mkdir /home/rstudio/Data
RUN mkdir /home/rstudio/Results

# ------------------
# --- Annotation ---
# ------------------
# Affymetrix
ADD Annotation/GPL570.txt /home/rstudio/Annotation/
ADD Annotation/GPL571.txt /home/rstudio/Annotation/
ADD Annotation/GPL96.txt /home/rstudio/Annotation/
ADD Annotation/GPL6244.txt /home/rstudio/Annotation/
ADD Annotation/GPL3535.txt /home/rstudio/Annotation/
ADD Annotation/GPL6246.txt /home/rstudio/Annotation/
ADD Annotation/GPL1261.txt /home/rstudio/Annotation/
ADD Annotation/GPL570.txt /home/rstudio/Annotation/
ADD Annotation/GPL1355.txt /home/rstudio/Annotation/

# ----------------------
# --- Pre Processing ---
# ----------------------
ADD Parameters.R /home/rstudio/

ADD Module-A/ModuleProcessingShiny.R /home/rstudio/Module-A/
ADD Module-A/Dependences.R /home/rstudio/Module-A/
ADD Module-A/ModuleProcessing.R /home/rstudio/Module-A/
ADD Module-A/cNodeExperiment.R /home/rstudio/Module-A/
ADD Module-A/cNodeValues.R /home/rstudio/Module-A/
ADD Module-A/Pre_Raw.R /home/rstudio/Module-A/
ADD Module-A/Pre_eSet.R /home/rstudio/Module-A/
ADD Module-A/Pre_Limma.R /home/rstudio/Module-A/
ADD Module-A/Pre_PosLimma.R /home/rstudio/Module-A/
ADD Module-A/MyStandardise.R /home/rstudio/Module-A/
ADD Module-A/Enrichment_GOStat2.R /home/rstudio/Module-A/
ADD Module-A/HeatMap3.R /home/rstudio/Module-A/
ADD Module-A/Pop_enrichment.R /home/rstudio/Module-A/
ADD Module-A/Clust_pop.R /home/rstudio/Module-A/
ADD Module-A/Clusterization.R /home/rstudio/Module-A/
ADD Module-A/GraphDataPersist.R /home/rstudio/Module-A/
ADD Module-A/PowerDetec.R /home/rstudio/Module-A/
ADD Module-A/AmdDetec.R /home/rstudio/Module-A/

USER root
CMD chmod -R 777 /home/rstudio && /usr/lib/rstudio-server/bin/rserver --server-daemonize=1 && \
    /usr/local/neo4j/bin/neo4j start && sleep 20 && /usr/bin/shiny-server.sh
