#!/usr/bin/bash -e

if [[ $EB_IS_COMMAND_LEADER == "true" ]];
then
	source "$PYTHONPATH/activate" && {
	# log which migrations have already been applied
	python manage.py showmigrations;
	# migrate
	python manage.py migrate --noinput;
	# bootstrap superuser
	python manage.py mysuperuser;
	# collect static 
	python manage.py collectstatic --noinput;
	}
else
  echo "this instance is NOT the leader";
fi
