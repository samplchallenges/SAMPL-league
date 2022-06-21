#!/bin/bash -e

sudo certbot -n -d app-staging.samplchallenges.org --nginx --agree-tos --email mhenry5@uci.edu
