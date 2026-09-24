"""Per-token prices used to calculate cost, in US$ per 1 million tokens: (input, output).

This is the ONE place to update prices. Output prices include thinking/reasoning tokens.
OpenRouter models aren't listed here: OpenRouter reports the cost of each request itself.
"""

PRICES_PER_MILLION = {
    # Anthropic API, standard pricing. Checked 2026-09-24.
    "claude-sonnet-5": (2.00, 10.00),
    "claude-opus-5": (5.00, 25.00),
    "claude-opus-5-5": (4.00, 20.00),
    "claude-haiku-4-5": (1.00, 5.00),
    # Gemini API (Google AI Studio), paid Standard tier, from ai.google.dev/gemini-api/docs/pricing.
    # Checked 2026-09-24. NOTE: rises to (1.50, 7.50) on 2027-01-01 - update this line then.
    "gemini-3.8-flash": (0.75, 3.75),
}


def cost_usd(model, input_tokens, output_tokens):
    """Cost of one request; NaN if the model has no price listed (so it's never silently $0)."""
    price_in, price_out = PRICES_PER_MILLION.get(model, (float("nan"), float("nan")))
    return (input_tokens * price_in + output_tokens * price_out) / 1e6
