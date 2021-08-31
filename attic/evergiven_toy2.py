import time


if __name__ == "__main__":
    import django
    django.setup()

from ever_given import wrapper
from core import models

def main():
    container_uri = "robbason/logging-example:latest"

    evaluation = models.Evaluation.objects.get(id=85)


    log_handler=models.Evaluation.LogHandler(evaluation)

    # final_command = wrapper._prepare_commandline("", {"smiles": "c1ccccc1"})
    results = wrapper.run(
        container_uri,
        file_kwargs={},
        kwargs={"smiles": "c1ccccc1"},
        log_handler=log_handler,
    )

    for k, v in results:
        print(f"Result: [{k}]", v)


    print("ALL DONE")

if __name__ == "__main__":
    main()
