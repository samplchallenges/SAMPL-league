import docker

client = docker.from_env()

BATCH = True

if BATCH:
    answers = []
    errors = []
    test_mols = [_ * "C" for _ in range(1, 10)]
    for test_mol in test_mols:
        result = client.containers.run(
            "mmh42/sampl-test:0.1", test_mol, auto_remove=True
        )
        result = float(result.strip())
        answers.append(result)
        print(test_mol, result)
    for answer, test_mol in zip(answers, test_mols):
        # I don't think options are being passed correctly here, all the values seem the same
        prediction = client.containers.run(
            "mmh42/sampl-test:0.1", test_mol, "--fuzz", auto_remove=True
        )
        prediction = float(prediction.strip())
        print(test_mol, prediction)
        errors.append(prediction - result)

    # RMSE
    RMSE = sum([_ ** 2 for _ in errors]) / len(errors)
    print(f"score: {RMSE}")

else:
    result = client.containers.run("mmh42/sampl-test:0.1", "CCCC", auto_remove=True)
    result = float(result.strip())
    print(result)
