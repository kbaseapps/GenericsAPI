FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN echo "start building docker image"

# R related installations
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
RUN echo 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list

RUN apt-get update --fix-missing
RUN apt-get install -y gcc wget r-base r-base-dev

RUN cp /usr/bin/R /kb/deployment/bin/.
RUN cp /usr/bin/Rscript /kb/deployment/bin/.

## Install packages are available for ecologists
# vegan: Community Ecology Package
RUN Rscript -e "install.packages('vegan')"

RUN pip install --upgrade pip \
    && python --version

RUN pip install numpy==1.19.5 \
    && pip install networkx==2.5.1 \
    && pip install pandas==1.1.5 \
    && pip install xlrd==2.0.1 \
    && pip install openpyxl==3.0.9 \
    && pip install xlsxwriter==3.0.2 \
    && pip install dotmap==1.3.25 \
    && pip install matplotlib==3.3.4 \
    && pip install scipy==1.5.4 \
    && pip install natsort==8.0.0 \
    && pip install scikit-learn==0.24.2 \
    && pip install plotly==5.3.1 \
    && pip install mock==4.0.3 \
    && pip install biom-format==2.1.10 \
    && pip install datashader==0.13.0 \
    && pip install rpy2==3.4.5 \
    && pip install coverage==6.1.1 \
    && pip install scikit-bio==0.5.6

# -----------------------------------------

RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

# install SigmaJS http://sigmajs.org/
RUN cd /kb/module && \
    wget https://github.com/jacomyal/sigma.js/archive/v1.2.1.zip && \
    unzip v1.2.1.zip && \
    rm -rf v1.2.1.zip && \
    mv sigma.js-1.2.1 sigma_js

# -----------------------------------------

ENV PYTHONUNBUFFERED=1

COPY ./ /kb/module
WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
