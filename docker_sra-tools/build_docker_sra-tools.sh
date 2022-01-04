cd /Users/robertshear/Projects/netseq/docker_sra-tools
imagename=rdshear/sra-tools
tagname=${imagename}:$(cat ../VERSION)
docker build . -t ${tagname}
docker tag ${tagname} ${imagename}:latest
docker push ${tagname}
docker push ${imagename}:latest
docker system prune -f