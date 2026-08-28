"""Help text for a form option: the sentence and links behind its (?) button.

Lives in the shared library, keyed by option id, so one entry serves the option
in every genome that offers it — eighteen options are declared by both
human_grch38.json and human_grch37.json.

An option offered at a different version per assembly — gnomAD is v4.1 on
GRCh38 and v2.1 on GRCh37 — is authored once: `{version}` stands in for the
version, and a link states the `major_version` it describes. Both are resolved
against the option's own label when the payload is built, so what is served is
finished text and the links that belong with it.
"""

import re

from pydantic import BaseModel, ConfigDict, Field

# A version in an option's label: "gnomAD Exomes v4.1.1" -> "v4.1.1".
_VERSION_IN_LABEL = re.compile(r"\bv\d+(?:\.\d+)*")
# Takes the preceding space with it, so a description reads cleanly when the
# label carries no version.
_VERSION_PLACEHOLDER = re.compile(r"\s?\{version\}")


class OptionHelpLink(BaseModel):
    """A resource link rendered after the description."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    href: str
    # Visible link text. The frontend supplies a generic label when unset.
    label: str | None = None
    # The major version this link describes. A link that states one is served
    # only to an option whose label carries that version, so an assembly cites
    # its own release rather than the other one's; a link without one is always
    # served. Consumed when the payload is built, never sent.
    major_version: str | None = Field(default=None, alias="majorVersion")


class OptionHelp(BaseModel):
    """One option's help, keyed by the id the form and results both use."""

    model_config = ConfigDict(extra="forbid")

    option_id: str
    # A `*span*` renders emphasised — a small markdown subset, so the text stays
    # a plain string rather than markup the backend has to escape.
    description: str
    links: list[OptionHelpLink] = []
    # Help that belongs to the input form alone. The results view reads the same
    # option payload and hangs an option's help on whichever node turns out to
    # be its visible title — which, for an option whose rows carry help of their
    # own, is a row that already has a (?) button. Set this where the rows say
    # everything the option would.
    form_only: bool = False

    def as_payload(self, label: str) -> dict:
        """The finished `help` object served on a form option.

        `label` is the option's own, and is what `{version}` and a link's
        `major_version` resolve against.

        `option_id` is the key rather than part of the value, and
        `major_version` has been spent choosing the links, so neither is sent;
        unset fields are dropped too, keeping a one-key link one key wide rather
        than padding every one of them with nulls.
        """
        version_match = _VERSION_IN_LABEL.search(label)
        version = version_match.group(0) if version_match else None
        # "v4.1" -> "4". Matching on the major alone keeps a link that still
        # describes the right callset through a point release.
        major = version[1:].split(".")[0] if version else None

        payload = self.model_dump(
            mode="json",
            by_alias=True,
            exclude_none=True,
            # `form_only` says where the help goes, not what it says, and is
            # acted on before this — never sent.
            exclude={"option_id", "links", "form_only"},
        )
        payload["description"] = _VERSION_PLACEHOLDER.sub(
            f" {version}" if version else "", self.description
        )
        links = [
            link.model_dump(
                mode="json", by_alias=True, exclude_none=True,
                exclude={"major_version"},
            )
            # A version-specific link is dropped rather than guessed at when the
            # label carries no version: citing the wrong release is worse than
            # citing none.
            for link in self.links
            if link.major_version is None or link.major_version == major
        ]
        if links:
            payload["links"] = links
        return payload


class HelpSpec(BaseModel):
    """The help half of the shared library."""

    model_config = ConfigDict(extra="forbid")

    options: list[OptionHelp] = []

    def form_only_options(self) -> set[str]:
        """Options whose help the results view should not show."""
        return {o.option_id for o in self.options if o.form_only}

    def payload_for(self, option_id: str, label: str) -> dict | None:
        """This option's finished help payload, or None when it has none.

        An option without help is ordinary — nothing is offered for
        `updownstream_distance` or the ClinVar sub-entries — so a miss is not an
        error and the key is simply omitted from the served option.
        """
        for option in self.options:
            if option.option_id == option_id:
                return option.as_payload(label)
        return None
