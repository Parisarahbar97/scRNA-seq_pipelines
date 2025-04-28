docker run -it -p 8888:8888 \
-v /var/lib/docker/parisa_tmp:/data \
ah3918/sc_analysis \
jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root