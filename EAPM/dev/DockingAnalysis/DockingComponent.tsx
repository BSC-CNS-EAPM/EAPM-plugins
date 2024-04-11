import { useState } from "react";

export default function DockingAnalysisView() {
  const [count, setCount] = useState(0);

  return (
    <div>
      <p>Hello World</p>
      <p>Test</p>
      <div>Counter: {count}</div>
      <button
        onClick={() => {
          setCount(count + 1);
        }}
      >
        Increase count
      </button>
      <div>Another div</div>
    </div>
  );
}
