FROM ubuntu:22.04

# Create a non-root user
RUN groupadd -r myuser && useradd -r -g myuser myuser

WORKDIR /app

RUN apt-get update
RUN apt install build-essential -y

RUN apt install gfortran -y

RUN apt install libnetcdf-dev -y
RUN apt install libnetcdff-dev -y

RUN apt install libopenmpi-dev -y

RUN apt-get clean

COPY . .

RUN nc-config --all
RUN nf-config --all

RUN mkdir /app/src/include

RUN cd /app/src && make

# Change the ownership of the application files to the non-root user
RUN chown -R myuser:myuser /app
# Switch to the non-root user
USER myuser

#RUN cd /app/src && ./main.exe ../config/namelist.in 12
RUN cd /app/src && mpiexec -n 8 ./main.exe ../config/namelist.in

CMD [ "sleep", "infinity" ]






