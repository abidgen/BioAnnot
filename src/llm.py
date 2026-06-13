"""Shared Anthropic tool-use call helper.

Single home for the forced-tool-use + prompt-caching call pattern that the
extractor and merger both repeat: build a cacheable system block, mark the tool
schema ephemeral, force ``tool_choice``, log token/cache usage, and return the
tool input. Importing this module never requires an API key (the client is built
lazily inside the call).
"""

from __future__ import annotations

import logging

import anthropic

log = logging.getLogger(__name__)


def get_client() -> anthropic.Anthropic:
    """Construct an Anthropic client (resolves ANTHROPIC_API_KEY from env)."""
    return anthropic.Anthropic()


async def call_tool(
    model: str,
    system_prompt: str,
    messages: list[dict],
    tools: list[dict],
    tool_name: str,
    max_tokens: int = 4096,
    label: str = "",
) -> tuple[dict, dict]:
    """Call Claude with forced tool use and prompt caching.

    The system prompt and every tool schema are marked ``cache_control:
    ephemeral`` so the static prefix is served from the prompt cache across
    calls. ``tool_choice`` forces ``tool_name`` so the model must emit exactly
    one structured tool call. Returns ``(tool_input_dict, usage_dict)``.
    """
    client = get_client()
    response = client.messages.create(
        model=model,
        max_tokens=max_tokens,
        system=[{
            "type": "text",
            "text": system_prompt,
            "cache_control": {"type": "ephemeral"}
        }],
        tools=[{**t, "cache_control": {"type": "ephemeral"}}
               for t in tools],
        tool_choice={"type": "tool", "name": tool_name},
        messages=messages,
    )
    usage = response.usage
    log.info(
        "%s tokens: input=%s output=%s cache_read=%s cache_created=%s",
        label,
        usage.input_tokens,
        usage.output_tokens,
        getattr(usage, "cache_read_input_tokens", 0),
        getattr(usage, "cache_creation_input_tokens", 0),
    )
    tool_block = next(
        (b for b in response.content if b.type == "tool_use"),
        None
    )
    if tool_block is None:
        raise ValueError(
            f"No tool_use block in response for {label}. "
            f"Content types: {[b.type for b in response.content]}"
        )
    return tool_block.input, usage.__dict__
