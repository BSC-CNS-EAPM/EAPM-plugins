import { createRoot } from "react-dom/client";
import { NBDSuiteResults } from "./nbdsuite_gui";

let container = null;

document.addEventListener("DOMContentLoaded", () => {
    if (!container) {
        container = document.getElementById("nbdsuite_root");
        const root = createRoot(container);
        root.render(<NBDSuiteResults />);
    }
});

