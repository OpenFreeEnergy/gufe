{% block title -%}

.. raw:: html

   <div style="display: None;">

{{ ("``" ~ objname ~ "``") | underline('=')}}

.. raw:: html

   </div>

{%- endblock %}
{% block base %}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :member-order: alphabetical  {# For consistency with Autosummary #}
   {% if show_inherited_members %}:inherited-members:
   {% endif %}{% if show_undoc_members %}:undoc-members:
   {% endif %}{% if show_inheritance %}:show-inheritance:
   {% endif %}

   {% block methods %}

   {%- set doc_methods = [] -%}
   {%- for item in methods -%}
      {%- if item not in ["__new__", "__init__"] and (show_inherited_members or item not in inherited_members) -%}
         {%- set _ = doc_methods.append(item) -%}
      {%- endif -%}
   {%- endfor %}

   {% if doc_methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :nosignatures:
   {% for item in doc_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}

   {%- set doc_attributes = [] -%}
   {%- for item in attributes -%}
      {%- if show_inherited_members or item not in inherited_members -%}
         {%- set _ = doc_attributes.append(item) -%}
      {%- endif -%}
   {%- endfor %}

   {% if doc_attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
      :nosignatures:
   {% for item in doc_attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
{% endblock %}
