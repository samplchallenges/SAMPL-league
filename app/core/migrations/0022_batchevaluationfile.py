# Generated by Django 4.0.6 on 2022-08-03 10:02

import core.models.run_related
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0021_alter_valuetype_batch_method'),
    ]

    operations = [
        migrations.CreateModel(
            name='BatchEvaluationFile',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('data', models.FileField(upload_to=core.models.run_related.batch_evaluation_upload_location)),
                ('batch_evaluation', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='core.batchevaluation')),
                ('value_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='core.valuetype')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]