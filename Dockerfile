FROM python:3.9

RUN apt-get install python3-tk

WORKDIR /code

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . .
