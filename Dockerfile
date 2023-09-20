FROM python:3.9

RUN apt-get install python-tk

WORKDIR /code

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY /src .

CMD ["python", "main.py"]