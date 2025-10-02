FILE=bfbrar

all: docker
.PHONY: all dbuild docker

## build Docker image (requires root access for Docker)
dbuild: Dockerfile
	docker build \
    -t $(FILE) .

## run RStudio in Docker container
docker-rstudio: dbuild
	echo "open RStudio Server at http://localhost:8787"
	docker run \
    --rm \
	-e DISABLE_AUTH=true \
	-e ROOT=true \
	-e USERID=$(id -u) \
	-e GROUPID=$(id -g) \
	-p 8787:8787 \
	-v $(CURDIR)/paper:/home/rstudio/paper \
	$(FILE)
##  --rm : automatically remove container when exiting
##  -e DISABLE_AUTH=true : disable RStudio Server authentication
##  -e ROOT=true : login as root user
##  -e USERID=$(id -u) : use host's user ID to match permissions
##  -e GROUPID=$(id -g) : use host's group ID to match permissions
##  -p 8787:8787 : map port 8787 of container to port 8787 of host
##  -v /$(CURDIR):/home/rstudio/paper : mount current directory into container
##  $(FILE) : run the $(FILE) image

## generate tex from Docker container and compile outside to PDF
docker-paper: dbuild
	docker run \
	--rm \
	--user rstudio \
	-e USERID=$(id -u) \
	-e GROUPID=$(id -g) \
	-v $(CURDIR)/paper:/home/rstudio/paper \
	-w /home/rstudio/paper \
	$(FILE) \
	make tex
	cd paper && make pdf2
