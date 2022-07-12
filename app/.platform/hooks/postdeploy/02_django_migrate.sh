#!/usr/bin/env bash

if [[ $EB_IS_COMMAND_LEADER == "true" ]];
then
	source "$PYTHONPATH/activate" && {
	# log which migrations have already been applied
	python manage.py showmigrations;
	# migrate
	python manage.py migrate --noinput;
	# bootstrap superuser
	python manage.py mysuperuser;
	# fill site with some toy data
	python manage.py sample_data --delete;
	# collect static
	python manage.py collectstatic --noinput;
	}
else
  echo "this instance is NOT the leader";
fi
