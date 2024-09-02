# https://hub.docker.com/r/library/python/tags/
#FROM python:3-onbuild  
FROM python:3.11.0

WORKDIR /usr/src/app/
# Copying source in current directory into the image
# python:3-onbuild expects the source in /usr/src/app
COPY biosimdb-app /usr/src/app/
RUN pip install -e .

# Run app with gunicorn
CMD ["gunicorn", "-c", "/usr/src/app/docker/gunicorn_config.py", "--bind", "0.0.0.0:80", "--forwarded-allow-ips='*'", "--chdir", "/usr/src/app/biosimdb_app", "app:app"] 

# Use line below for development mode
# CMD ["flask", "--app", "biosimdb_app", "run", "--debug"]
