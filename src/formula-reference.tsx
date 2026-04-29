import { Action, ActionPanel, Detail, showToast, Toast } from "@raycast/api";
import { useEffect, useState } from "react";
import { MATH_101 } from "./math/content";
import { clearMathCache } from "./math/render";
import { processMath } from "./math/process";

export default function FormulaReference() {
  const [markdown, setMarkdown] = useState<string>("");
  const [loading, setLoading] = useState(true);
  const [reloadKey, setReloadKey] = useState(0);

  useEffect(() => {
    let cancelled = false;
    (async () => {
      setLoading(true);
      try {
        const rendered = await processMath(MATH_101);
        if (!cancelled) setMarkdown(rendered);
      } catch (e) {
        if (!cancelled)
          setMarkdown(
            `# Error\n\n\`\`\`\n${e instanceof Error ? e.message : String(e)}\n\`\`\``,
          );
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => {
      cancelled = true;
    };
  }, [reloadKey]);

  return (
    <Detail
      markdown={markdown}
      isLoading={loading}
      navigationTitle="Formula Reference"
      actions={
        <ActionPanel>
          <Action
            title="Reload"
            onAction={() => setReloadKey((k) => k + 1)}
            shortcut={{ modifiers: ["cmd"], key: "r" }}
          />
          <Action
            title="Clear Math Cache & Reload"
            onAction={async () => {
              const n = await clearMathCache();
              await showToast({
                style: Toast.Style.Success,
                title: `Cleared ${n} cached formulas`,
              });
              setReloadKey((k) => k + 1);
            }}
            shortcut={{ modifiers: ["cmd", "shift"], key: "r" }}
          />
        </ActionPanel>
      }
    />
  );
}
