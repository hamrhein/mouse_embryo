FROM debian:buster
MAINTAINER Henry Amrhein (hamrhein@caltech.edu)

ENV PYTHONPATH="/software/sdbl"
ENV PATH="/software/hamrhein:${PATH}"

RUN echo "deb http://ftp.us.debian.org/debian/ buster contrib non-free" >> /etc/apt/sources.list.d/01_repos.list

RUN apt-get update && apt-get upgrade -y

RUN apt-get install -y apt-utils

RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN echo ttf-mscorefonts-installer msttcorefonts/present-mscorefonts-eula note | debconf-set-selections

RUN apt-get install -y unzip wget python3 python3-numpy python3-pandas python3-matplotlib python3-pygraphviz python3-imageio ttf-mscorefonts-installer
RUN apt-get clean

COPY *.py /software/hamrhein/
COPY util/*.py util/*.sh /software/hamrhein/

WORKDIR /data

RUN wget -O MouseLimbData.h5 https://woldlab.caltech.edu/nextcloud/index.php/s/syNtQbdGessF5NB/download
RUN wget -O peng_bloom.zip https://woldlab.caltech.edu/nextcloud/index.php/s/MQ7DWssYTfmmPnR/download
RUN unzip peng_bloom.zip && rm peng_bloom.zip

WORKDIR /software/hamrhein
RUN chmod 0755 *

WORKDIR /output
CMD ["/software/hamrhein/runscript.sh"]
