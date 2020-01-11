FROM debian:buster
MAINTAINER Henry Amrhein (hamrhein@caltech.edu)

ENV PYTHONPATH="/software/sdbl"
ENV PATH="/software/hamrhein:${PATH}"

RUN echo "deb http://ftp.us.debian.org/debian/ buster contrib non-free" >> /etc/apt/sources.list.d/01_repos.list

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get install -y apt-utils

RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN echo ttf-mscorefonts-installer msttcorefonts/present-mscorefonts-eula note | debconf-set-selections

RUN apt-get install -y unzip wget python3 python3-numpy python3-pandas python3-matplotlib python3-pygraphviz python3-imageio ttf-mscorefonts-installer
RUN apt-get clean

WORKDIR /data

RUN wget -O MouseLimbData.h5 https://woldlab.caltech.edu/nextcloud/index.php/s/3nMkLMckn7Gtzdr/download
RUN wget -O peng_bloom.zip https://woldlab.caltech.edu/nextcloud/index.php/s/8kZ7dPXrMPXnEAJ/download
RUN unzip peng_bloom.zip
RUN rm peng_bloom.zip

WORKDIR /software/sdbl

RUN wget https://github.com/hamrhein/sdbl/archive/master.zip
RUN unzip -j master.zip
RUN rm master.zip

RUN mkdir /software/hamrhein
RUN mv build_sql_stringdb_database.py /software/hamrhein
RUN mv build_10x_tf_graphs.py /software/hamrhein
RUN mv build_blossom_graph.py /software/hamrhein
RUN mv runscript.sh /software/hamrhein

WORKDIR /software/hamrhein
RUN chmod 0755 *.py
RUN chmod 0755 runscript.sh

WORKDIR /output
CMD ["/software/hamrhein/runscript.sh"]
