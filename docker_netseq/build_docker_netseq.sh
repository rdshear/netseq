cd /Users/robertshear/Projects/netseq/docker_netseq
imagename=rdshear/netseq
tagname=${imagename}:$(cat ../VERSION)
docker build . -t ${tagname}
docker push ${tagname}