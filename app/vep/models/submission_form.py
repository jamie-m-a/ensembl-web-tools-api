from pydantic import BaseModel

class DropdownOption(BaseModel):
  label: str
  value: str

class Dropdown(BaseModel):
    label: str
    description: str | None = None
    type: str = "select"
    options: list[DropdownOption]
    default_value: str

class FormConfig(BaseModel):
    transcript_set: Dropdown

class GenomeAnnotationProvider(BaseModel):
  annotation_provider_name: str
  annotation_version: str
