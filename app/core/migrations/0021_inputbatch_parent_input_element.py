# Generated by Django 4.0.5 on 2022-06-16 01:05

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0020_filevalue_input_element_floatvalue_input_element_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='inputbatch',
            name='parent_input_element',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='core.inputelement'),
        ),
    ]
