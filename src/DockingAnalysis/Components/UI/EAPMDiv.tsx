import { HTMLAttributes } from "react";

export default function EAPMDiv(props: HTMLAttributes<HTMLDivElement>) {
  return (
    <div
      {...props}
      className={`${props.className} border-blue-500 border-2 shadow-lg shadow-blue-500/50 dark:shadow-lg font-medium rounded-lg text-sm px-5 py-2.5 text-center me-2 mb-2`}
    >
      {props.children}
    </div>
  );
}
