# the gitlab testing environment relies on docker image.
# If the test failed online and the error is hard to debug, it
# might be useful to setup a local docker container to run the 
# test within in order to reproduce the error. In this case, follow
# the instruction below.


# pull docker image with python 3.7 (debian stretch)
docker pull python:3.8-buster

# run bash shell on the docker image with a shared folder
docker run -it -v /media/jkl/data/phd/tmp/resipy:/resipy python:3.8-buster bash

# this will create a docker image in which you will need to run
# the commands in .gitlab-ci.yml in order to install necessary
# packages, then you can run the test `python test.py`




# for BINDER testing, we need the latest image of ubuntu:
docker pull ubuntu:latest

# and install
sudo apt-get install python3 pip3


# run existing container
docker start <container name or id>
docker exec -it <container name or id> bash

# run jupyterlab inside container
python -m jupyterlab --allow-root --ip 0.0.0.0



