"""Claude models go straight to the Anthropic API (ANTHROPIC_API_KEY), not OpenRouter."""
import anthropic

# $ per million tokens (input, output). Thinking tokens are billed as output.
PRICES = {
    "claude-sonnet-5": (2.00, 10.00),
    "claude-opus-5": (5.00, 25.00),
    "claude-opus-5-5": (4.00, 20.00),
    "claude-haiku-4-5": (1.00, 5.00),
}

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
    price_in, price_out = PRICES.get(model, (float("nan"), float("nan")))
    usage = {
        "prompt_tokens": resp.usage.input_tokens,
        "completion_tokens": resp.usage.output_tokens,  # includes thinking
        "reasoning_tokens": None,  # the API doesn't split thinking out of output_tokens
        "cost_usd": (resp.usage.input_tokens * price_in + resp.usage.output_tokens * price_out) / 1e6,
    }
    return text, resp.stop_reason, usage
