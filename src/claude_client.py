"""Claude models go straight to the Anthropic API (ANTHROPIC_API_KEY), not OpenRouter."""
import anthropic

from prices import cost_usd

client = anthropic.Anthropic(timeout=120, max_retries=6)


def ask_claude(model, prompt, max_tokens):
    """Return (text, stop_reason, usage dict). Sonnet 5 / Opus 5 reject `temperature`,
    and think adaptively by default, so neither is set here."""
    resp = client.messages.create(
        model=model,
        max_tokens=max_tokens,
        messages=[{"role": "user", "content": prompt}],
    )
    text = "".join(b.text for b in resp.content if b.type == "text")
    usage = {
        "prompt_tokens": resp.usage.input_tokens,
        "completion_tokens": resp.usage.output_tokens,  # includes thinking
        "reasoning_tokens": None,  # the API doesn't split thinking out of output_tokens
        "cost_usd": cost_usd(model, resp.usage.input_tokens, resp.usage.output_tokens),
    }
    return text, resp.stop_reason, usage
