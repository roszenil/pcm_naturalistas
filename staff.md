---
layout: page
title: Instructores
description: Los instructores del Taller en Quito 2026
---

## Instructores

{% assign instructors = site.staffers | where: 'role', 'Instructora' %}
{% for staffer in instructors %}
{{ staffer }}
{% endfor %}

{% assign teaching_assistants = site.staffers | where: 'role', 'Teaching Assistant' %}
{% assign num_teaching_assistants = teaching_assistants | size %}
{% if num_teaching_assistants != 0 %}

{% for staffer in teaching_assistants %}
{{ staffer }}
{% endfor %}
{% endif %}
