cd /Users/robertshear/Projects/netseq/docker_bbtools
imagename=rdshear/bbtools
tagname=${imagename}:$(cat ../VERSION)
docker build . --progress plain -t ${tagname}
docker tag ${tagname} ${imagename}:latest
docker push ${tagname}
docker push ${imagename}:latest
docker system prune -f