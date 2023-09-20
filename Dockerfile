FROM python:3.9

RUN sudo apt-get install python-tk

WORKDIR /code

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY . .

CMD ["python", "main.py"]
