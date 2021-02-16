from .settings import *

DEBUG = False
SECRET_KEY = os.environ["SAMPL_SECRET_KEY"]
ALLOWED_HOSTS = [
    "app.samplchallenges.org",
    "samplmvp-env.eba-rhcwa63p.us-east-2.elasticbeanstalk.com",
]
