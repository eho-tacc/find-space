FROM sd2e/python3:ubuntu18
LABEL maintainer="eho@tacc.utexas.edu"

### DO NOT EDIT BELOW ###
# Install software to make containers interactive
#   Debian-like
ARG FOUNDATIONS="bash curl git man rsync vim"
RUN if [ -f "/etc/os-release" ] && \
    grep "debian" /etc/os-release; \
    then apt-get update -qqq && \
    apt-get install -qqq -y -u apt-utils ${FOUNDATIONS} && \
    apt-get clean -qqq; fi
#   RHEL-like
RUN if [ -f "/etc/os-release" ] && \
    grep "rhel" /etc/os-release; then \
    yum -q updateinfo && \
    yum -q -y install ${FOUNDATIONS} && \
    yum -q clean all; rm -rf /var/cache/yum; fi
#   Alpine
RUN if [ -f "/etc/os-release" ] && \
    grep "alpine" /etc/os-release; then \
    apk update -q && \
    apk add -q ${FOUNDATIONS}; \
    fi
# Add shadow mount points image mapped to TACC storage
RUN for MOUNT in corral data gpfs projects scratch work; \
    do mkdir -p /${MOUNT}; \
    chown root:root /${MOUNT}; \
    chmod ug+rwx,o+rx /${MOUNT}; \
    done
# These variables used by TACC.cloud components
ENV _PROJ_CORRAL=/corral
ENV _PROJ_STOCKYARD=/work/projects/SD2E-Community
ENV _USER_WORK=
### DO NOT EDIT ABOVE ###

# Customize your Docker container starting here

USER root
RUN apt-get -y update && \
    apt-get -y install python3.6-dev
RUN pip3 install numpy \
                 scipy \
                 cython

# install ProtPy
RUN pip3 install --index-url https://test.pypi.org/simple/ protpy

ARG SRC="/opt/find-space/src"
ARG WD="/home/work"
RUN mkdir -p $SRC $WD && \
    cd $WD && \
    mkdir -p out scan log && \
    echo 1crn > pdb_screen_list.txt
COPY src $SRC

#CMD "bash"
CMD python3 /opt/find-space/src/find_space.py -w /home/work -p /home/work/pdb_screen_list.txt -l /home/work/log/fs_log.txt
