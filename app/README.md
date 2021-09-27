# This is a work in progress and is experimental

Demo site [here](http://samplmvp-env.eba-rhcwa63p.us-east-2.elasticbeanstalk.com)

## Getting started locally

### Postgres
- install postgres locally
- create a postgres database, say 'sampl'

### Pipenv
`install pipenv`
`pipenv install`
`pipenv shell` (activates your python environment)

### Environment variables
> `export RDS_DB_NAME=sampl`
> `export RDS_USERNAME=your_username`
> `export RDS_PASSWORD=your_password`
> `export DJANGO_SETTINGS_MODULE=sampl.settings_dev`

recommended to use direnv: https://direnv.net/ to help with these settings

### Django setup
`django-admin migrate` (will create necessary database tables)

If you want sample data (recommended to start):
`django-admin sample_data`
`django-admin createsuperuser`

### Running
Now, we need three terminals, each with the right environment variables and all with an active pipenv shell. In one:
`django-admin runserver`

In another:
`dask-scheduler`

In another:
`dask-worker --preload daskworkerinit.py  127.0.0.1:8786`

Then, in your browser, go to http://localhost:8080
