FROM ubuntu:14.04

RUN cp /etc/apt/sources.list /etc/apt/sources.list.bak
RUN echo " \n\
deb http://mirrors.aliyun.com/ubuntu/ trusty main restricted universe multiverse \n\
deb http://mirrors.aliyun.com/ubuntu/ trusty-security main restricted universe multiverse \n\
deb http://mirrors.aliyun.com/ubuntu/ trusty-updates main restricted universe multiverse \n\
deb http://mirrors.aliyun.com/ubuntu/ trusty-proposed main restricted universe multiverse \n\
deb http://mirrors.aliyun.com/ubuntu/ trusty-backports main restricted universe multiverse \n\
deb-src http://mirrors.aliyun.com/ubuntu/ trusty main restricted universe multiverse \n\
deb-src http://mirrors.aliyun.com/ubuntu/ trusty-security main restricted universe multiverse \n\
deb-src http://mirrors.aliyun.com/ubuntu/ trusty-updates main restricted universe multiverse \n\
deb-src http://mirrors.aliyun.com/ubuntu/ trusty-proposed main restricted universe multiverse \n\
deb-src http://mirrors.aliyun.com/ubuntu/ trusty-backports main restricted universe multiverse \n\
" >/etc/apt/sources.list

RUN apt-get update && apt-get install -y ca-certificates software-properties-common

# install python
RUN apt-get install -y python-pip python2.7 python2.7-dev

# install pypy
RUN add-apt-repository ppa:pypy/ppa
RUN apt-get update
RUN apt-get install -y pypy pypy-dev

# ADD after.py
RUN mkdir -p /tools/after
COPY . /tools/after
RUN chmod +x /tools/after/after.py
WORKDIR /tools/after
RUN make

ENTRYPOINT ["/usr/bin/pypy"]
